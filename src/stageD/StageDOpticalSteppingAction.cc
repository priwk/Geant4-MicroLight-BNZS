#include "StageDOpticalSteppingAction.hh"

#include "AnalysisConfig.hh"
#include "DetectorConstruction.hh"
#include "StageDOpticalEventAction.hh"
#include "StageDOpticalRunAction.hh"
#include "StageDReentrySampler.hh"

#include "G4EventManager.hh"
#include "G4Exception.hh"
#include "G4DynamicParticle.hh"
#include "G4OpticalPhoton.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4StackManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4VProcess.hh"
#include "G4VPhysicalVolume.hh"

#include <algorithm>
#include <cmath>
#include <string>

namespace
{
  constexpr G4double kBoundaryEpsilon = 1.0e-4 * um;

  G4bool IsOpticalPhoton(const G4Track *track)
  {
    return track != nullptr &&
           track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition();
  }

  G4double AngleDeg(const G4ThreeVector &a, const G4ThreeVector &b)
  {
    const G4double dot = std::clamp(a.dot(b), -1.0, 1.0);
    return std::acos(dot) / deg;
  }

  G4double ClampCosTheta(const G4ThreeVector &a, const G4ThreeVector &b)
  {
    return std::clamp(a.dot(b), -1.0, 1.0);
  }

  std::string ProcessName(const G4StepPoint *point)
  {
    if (point == nullptr)
      return "";
    const auto *process = point->GetProcessDefinedStep();
    if (process == nullptr)
      return "";
    return process->GetProcessName();
  }

  G4bool StartsWith(const G4String &s, const char *prefix)
  {
    const std::string value = s;
    const std::string p = prefix;
    return value.rfind(p, 0) == 0;
  }

  DetectorConstruction::Phase PhaseFromPhysicalVolume(const G4VPhysicalVolume *pv)
  {
    if (pv == nullptr)
      return DetectorConstruction::Phase::World;

    const auto *lv = pv->GetLogicalVolume();
    if (lv == nullptr)
      return DetectorConstruction::Phase::World;

    const auto &lvName = lv->GetName();
    if (lvName == "BN_LV" || StartsWith(lvName, "BN_ClipLV_"))
      return DetectorConstruction::Phase::BN;
    if (lvName == "ZnS_LV" || StartsWith(lvName, "ZnS_ClipLV_"))
      return DetectorConstruction::Phase::ZnS;
    if (lvName == "MatrixLV")
      return DetectorConstruction::Phase::Matrix;
    if (lvName == "WorldLV")
      return DetectorConstruction::Phase::World;
    return DetectorConstruction::Phase::Unknown;
  }

  G4bool IsOutsideRve(const G4ThreeVector &position,
                      const DetectorConstruction *detector)
  {
    if (detector == nullptr)
      return false;

    const G4double halfX = detector->GetPatchHalfXUm() * um;
    const G4double halfY = detector->GetPatchHalfYUm() * um;
    const G4double halfZ = detector->GetPatchHalfZUm() * um;
    return (std::abs(position.x()) > halfX ||
            std::abs(position.y()) > halfY ||
            std::abs(position.z()) > halfZ);
  }

  G4bool IsParticlePhase(DetectorConstruction::Phase phase)
  {
    return phase == DetectorConstruction::Phase::BN ||
           phase == DetectorConstruction::Phase::ZnS;
  }

  G4bool IsReentryPhase(DetectorConstruction::Phase phase)
  {
    return IsParticlePhase(phase) ||
           phase == DetectorConstruction::Phase::Matrix;
  }

  std::string PhaseLabel(DetectorConstruction::Phase phase)
  {
    if (phase == DetectorConstruction::Phase::Unknown)
      return "Unknown";
    return DetectorConstruction::PhaseName(phase);
  }

  G4bool ComputeExitPointOnRveBox(const G4ThreeVector &prePos,
                                  const G4ThreeVector &oldDir,
                                  const DetectorConstruction *detector,
                                  G4ThreeVector &exitPoint)
  {
    if (detector == nullptr || oldDir.mag2() <= 0.0)
      return false;

    const G4double halfX = detector->GetPatchHalfXUm() * um;
    const G4double halfY = detector->GetPatchHalfYUm() * um;
    const G4double halfZ = detector->GetPatchHalfZUm() * um;

    G4double tMin = -std::numeric_limits<G4double>::infinity();
    G4double tMax = std::numeric_limits<G4double>::infinity();

    auto updateAxis = [&](G4double pos, G4double dir, G4double half) -> G4bool
    {
      if (std::abs(dir) <= 1.0e-18)
        return std::abs(pos) <= half;

      G4double t1 = (-half - pos) / dir;
      G4double t2 = (+half - pos) / dir;
      if (t1 > t2)
        std::swap(t1, t2);
      tMin = std::max(tMin, t1);
      tMax = std::min(tMax, t2);
      return tMin <= tMax;
    };

    if (!updateAxis(prePos.x(), oldDir.x(), halfX) ||
        !updateAxis(prePos.y(), oldDir.y(), halfY) ||
        !updateAxis(prePos.z(), oldDir.z(), halfZ))
    {
      return false;
    }

    if (tMax <= 0.0 || !std::isfinite(tMax))
      return false;

    exitPoint = prePos + tMax * oldDir;
    return true;
  }

  StageDReentryDiagnosticRecord MakeReentryDiagnosticRecord(
      const StageDReentrySampler::ReentryContext &ctx,
      const StageDReentrySampler::ReentryDiagnostics &diag)
  {
    StageDReentryDiagnosticRecord record;
    record.event_id = ctx.eventID;
    record.reentry_index = ctx.reentryIndex;
    record.strategy = diag.strategy;
    record.fallback_level = diag.fallbackLevel;
    record.exit_phase = PhaseLabel(diag.exitPhase);
    record.entry_phase = PhaseLabel(diag.entryPhase);
    record.old_dir = ctx.oldDir;
    record.exit_point = ctx.exitPoint;
    record.entry_point = diag.entryPoint;
    record.particle_q_exit = diag.particleQExit;
    record.particle_q_entry = diag.particleQEntry;
    record.particle_mu_exit = diag.particleMuExit;
    record.particle_mu_entry = diag.particleMuEntry;
    record.matrix_clearance_exit_um = diag.matrixClearanceExitUm;
    record.matrix_clearance_entry_um = diag.matrixClearanceEntryUm;
    record.matrix_nearest_phase_exit = diag.matrixNearestPhaseExit;
    record.matrix_nearest_phase_entry = diag.matrixNearestPhaseEntry;
    record.matrix_clearance_bin_exit = diag.matrixClearanceBinExit;
    record.matrix_clearance_bin_entry = diag.matrixClearanceBinEntry;
    record.trials = diag.trials;
    return record;
  }

  void AccumulateReentryCounters(StageDPhotonEventRecord &event,
                                 const StageDReentrySampler::ReentryDiagnostics &diag)
  {
    if (diag.strategy == "particle_sphere_q_mu")
      ++event.num_reentry_particle_q_mu;
    else if (diag.strategy == "matrix_clearance_binned_portal")
      ++event.num_reentry_matrix_clearance_portal;
    else if (diag.strategy == "matrix_random_debug")
      ++event.num_reentry_random_matrix_debug;

    if (diag.fallbackLevel == "same_bin")
      ++event.num_reentry_fallback_same_bin;
    else if (diag.fallbackLevel == "adjacent_bin")
      ++event.num_reentry_fallback_adjacent_bin;
    else if (diag.fallbackLevel == "same_phase_any_bin")
      ++event.num_reentry_fallback_any_bin;
    else if (diag.fallbackLevel == "any_phase_same_bin")
      ++event.num_reentry_fallback_any_phase_same_bin;
    else if (diag.fallbackLevel == "any_portal")
      ++event.num_reentry_fallback_any_portal;
  }

  G4int PhaseFunctionBin(const G4double cosTheta)
  {
    constexpr G4double kMin = -1.0;
    constexpr G4double kMax = 1.0;
    constexpr G4double kSpan = kMax - kMin;
    const G4double shifted = std::clamp(cosTheta, kMin, kMax) - kMin;
    const G4double normalized = shifted / kSpan;
    const G4int index = static_cast<G4int>(
        normalized * static_cast<G4double>(StageDPhotonEventRecord::kPhaseFunctionBins));
    return std::clamp(index, 0,
                      static_cast<G4int>(StageDPhotonEventRecord::kPhaseFunctionBins) - 1);
  }
}

StageDOpticalSteppingAction::StageDOpticalSteppingAction(
    StageDOpticalRunAction *runAction,
    StageDOpticalEventAction *eventAction,
    AnalysisConfig *config)
    : G4UserSteppingAction(),
      fConfig(config),
      fRunAction(runAction),
      fEventAction(eventAction),
      fDetector(nullptr),
      fReentrySampler(nullptr)
{
}

StageDOpticalSteppingAction::~StageDOpticalSteppingAction()
{
  delete fReentrySampler;
}

const DetectorConstruction *StageDOpticalSteppingAction::ResolveDetector() const
{
  if (fDetector == nullptr)
  {
    fDetector = dynamic_cast<const DetectorConstruction *>(
        G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  }
  return fDetector;
}

G4bool StageDOpticalSteppingAction::HandleBoundaryReentry(
    const G4Step *step,
    G4Track *track,
    const DetectorConstruction *detector)
{
  if (fConfig == nullptr || fEventAction == nullptr || detector == nullptr)
    return false;

  auto &event = fEventAction->MutableCurrentEvent();
  const G4ThreeVector postPos = step->GetPostStepPoint()->GetPosition();
  const G4bool outOfBox = IsOutsideRve(postPos, detector);
  if (!outOfBox)
    return false;

  if (fConfig->stageD_boundary_mode == "escape")
  {
    event.in_particle_segment = false;
    event.particle_segment_phase.clear();
    fEventAction->SetFinalStatus("escaped_debug", false);
    track->SetTrackStatus(fStopAndKill);
    return true;
  }

  if (fConfig->stageD_boundary_mode != "same_phase_reentry")
    return false;

  if (event.num_reentry >= fConfig->stageD_max_reentry)
  {
    event.in_particle_segment = false;
    event.particle_segment_phase.clear();
    fEventAction->SetFinalStatus("max_reentry", false);
    track->SetTrackStatus(fStopAndKill);
    return true;
  }

  if (fReentrySampler == nullptr)
  {
    fReentrySampler = new StageDReentrySampler(detector, fConfig);
    if (fRunAction != nullptr)
      fRunAction->SetReentryPortalSummary(fReentrySampler->GetPortalSummary());
  }

  const G4ThreeVector prePos = step->GetPreStepPoint()->GetPosition();
  const G4ThreeVector oldDir = step->GetPreStepPoint()->GetMomentumDirection();
  G4ThreeVector exitPoint = prePos;
  ComputeExitPointOnRveBox(prePos, oldDir, detector, exitPoint);
  const G4ThreeVector exitInsidePoint = exitPoint - kBoundaryEpsilon * oldDir.unit();

  auto phase = fReentrySampler->FastPhaseAtPointForReentry(exitInsidePoint);
  if (phase == DetectorConstruction::Phase::World ||
      phase == DetectorConstruction::Phase::Unknown)
  {
    phase = fReentrySampler->FastPhaseAtPointForReentry(prePos);
  }

  StageDReentrySampler::ReentryContext ctx;
  ctx.phase = phase;
  ctx.prePos = prePos;
  ctx.postPos = postPos;
  ctx.oldDir = oldDir;
  ctx.exitPoint = exitPoint;
  ctx.exitInsidePoint = exitInsidePoint;
  ctx.wavelengthNm = event.wavelength_nm;
  ctx.eventID = event.photonID;
  ctx.reentryIndex = event.num_reentry + event.num_reentry_failed + 1;

  StageDReentrySampler::ReentryDiagnostics diag;
  diag.exitPhase = phase;
  diag.exitInsidePoint = exitInsidePoint;

  if (!IsReentryPhase(phase))
  {
    diag.strategy = "phase_resolve_failed";
    diag.fallbackLevel = "unsupported_exit_phase";
    if (fRunAction != nullptr)
      fRunAction->RecordReentryDiagnostic(MakeReentryDiagnosticRecord(ctx, diag));
    ++event.num_reentry_failed;
    event.in_particle_segment = false;
    event.particle_segment_phase.clear();
    fEventAction->SetFinalStatus("reentry_failed", false);
    track->SetTrackStatus(fStopAndKill);
    return true;
  }

  G4ThreeVector newPosition;
  const G4bool ok = fReentrySampler->SampleReentry(ctx, newPosition, diag);
  if (fRunAction != nullptr)
    fRunAction->RecordReentryDiagnostic(MakeReentryDiagnosticRecord(ctx, diag));

  if (!ok)
  {
    ++event.num_reentry_failed;
    event.in_particle_segment = false;
    event.particle_segment_phase.clear();
    fEventAction->SetFinalStatus("reentry_failed", false);
    track->SetTrackStatus(fStopAndKill);
    return true;
  }

  auto *dynamicParticle = new G4DynamicParticle(
      G4OpticalPhoton::OpticalPhotonDefinition(),
      oldDir,
      track->GetKineticEnergy());
  dynamicParticle->SetPolarization(track->GetPolarization());

  auto *continuationTrack = new G4Track(
      dynamicParticle,
      track->GetGlobalTime(),
      newPosition);
  continuationTrack->SetLocalTime(track->GetLocalTime());
  continuationTrack->SetProperTime(track->GetProperTime());
  continuationTrack->SetWeight(track->GetWeight());
  continuationTrack->SetParentID(track->GetParentID());

  auto *stackManager = G4EventManager::GetEventManager()->GetStackManager();
  if (stackManager == nullptr)
  {
    delete continuationTrack;
    ++event.num_reentry_failed;
    event.in_particle_segment = false;
    event.particle_segment_phase.clear();
    fEventAction->SetFinalStatus("reentry_failed", false);
    track->SetTrackStatus(fStopAndKill);
    return true;
  }
  stackManager->PushOneTrack(continuationTrack);

  ++event.num_reentry;
  if (phase == DetectorConstruction::Phase::BN)
    ++event.num_reentry_BN;
  else if (phase == DetectorConstruction::Phase::ZnS)
    ++event.num_reentry_ZnS;
  else if (phase == DetectorConstruction::Phase::Matrix)
    ++event.num_reentry_matrix;
  AccumulateReentryCounters(event, diag);

  event.in_particle_segment = false;
  event.particle_segment_phase.clear();
  track->SetTrackStatus(fStopAndKill);
  return true;
}

G4bool StageDOpticalSteppingAction::HandleLimitKills(const G4Step *step, G4Track *track)
{
  if (fConfig == nullptr || fEventAction == nullptr)
    return false;

  auto &event = fEventAction->MutableCurrentEvent();

  if (event.num_encounter_total >= fConfig->stageD_target_primary_scatter)
  {
    fEventAction->SetFinalStatus("target_primary_scatter", false);
    track->SetTrackStatus(fStopAndKill);
    return true;
  }

  if (event.num_steps >= fConfig->stageD_max_steps)
  {
    fEventAction->SetFinalStatus("max_steps", false);
    track->SetTrackStatus(fStopAndKill);
    return true;
  }

  if (event.total_path_length_um >= fConfig->stageD_max_path_length_um)
  {
    fEventAction->SetFinalStatus("max_path_length", false);
    track->SetTrackStatus(fStopAndKill);
    return true;
  }

  const std::string processName = ProcessName(step->GetPostStepPoint());
  if (processName == "OpAbsorption")
  {
    const auto *prePoint = step->GetPreStepPoint();
    const auto *prePV = (prePoint != nullptr) ? prePoint->GetPhysicalVolume() : nullptr;
    const auto phase = PhaseFromPhysicalVolume(prePV);
    fEventAction->MarkAbsorbed(DetectorConstruction::PhaseName(phase));
    return false;
  }

  if (track->GetTrackStatus() == fStopAndKill &&
      event.final_status == "in_progress")
  {
    fEventAction->SetFinalStatus("lost", false);
  }

  return false;
}

void StageDOpticalSteppingAction::UserSteppingAction(const G4Step *step)
{
  if (step == nullptr || fEventAction == nullptr)
    return;

  G4Track *track = step->GetTrack();
  if (!IsOpticalPhoton(track))
    return;

  auto &event = fEventAction->MutableCurrentEvent();
  ++event.num_steps;
  event.total_path_length_um += step->GetStepLength() / um;

  const auto *prePoint = step->GetPreStepPoint();
  const auto *postPoint = step->GetPostStepPoint();
  if (prePoint == nullptr || postPoint == nullptr)
    return;

  const auto *detector = ResolveDetector();
  const G4bool isOuterRveExit =
      (detector != nullptr && IsOutsideRve(postPoint->GetPosition(), detector));
  const auto *prePV = prePoint->GetPhysicalVolume();
  const auto *postPV = postPoint->GetPhysicalVolume();
  const auto prePhase = PhaseFromPhysicalVolume(prePV);
  const auto postPhase =
      isOuterRveExit ? DetectorConstruction::Phase::World
                     : PhaseFromPhysicalVolume(postPV);

  const G4double stepLengthUm = step->GetStepLength() / um;
  if (prePhase == DetectorConstruction::Phase::BN)
    event.path_length_bn_um += stepLengthUm;
  else if (prePhase == DetectorConstruction::Phase::ZnS)
    event.path_length_zns_um += stepLengthUm;
  else if (prePhase == DetectorConstruction::Phase::Matrix)
    event.path_length_matrix_um += stepLengthUm;
  else
    event.path_length_world_um += stepLengthUm;

  const std::string processName = ProcessName(postPoint);
  const G4bool isMaterialBoundary =
      !isOuterRveExit &&
      (processName == "OpBoundary" ||
       (prePhase != DetectorConstruction::Phase::Unknown &&
        postPhase != DetectorConstruction::Phase::Unknown &&
        prePhase != postPhase &&
        prePhase != DetectorConstruction::Phase::World &&
        postPhase != DetectorConstruction::Phase::World));
  if (isMaterialBoundary)
  {
    ++event.num_material_boundary;
  }

  const G4ThreeVector preDir = prePoint->GetMomentumDirection();
  const G4ThreeVector postDir = postPoint->GetMomentumDirection();

  if (IsParticlePhase(prePhase) &&
      (!event.in_particle_segment || event.particle_segment_phase != DetectorConstruction::PhaseName(prePhase)))
  {
    event.in_particle_segment = true;
    event.particle_segment_phase = DetectorConstruction::PhaseName(prePhase);
    event.particle_segment_entry_direction = preDir;
  }

  const G4bool isParticleExitInsideRve =
      event.in_particle_segment &&
      IsParticlePhase(prePhase) &&
      prePhase != postPhase &&
      !isOuterRveExit &&
      postPhase != DetectorConstruction::Phase::World;

  if (isParticleExitInsideRve)
  {
    const G4double cosTheta = ClampCosTheta(
        event.particle_segment_entry_direction,
        postDir);
    const G4double oneMinusCosTheta = 1.0 - cosTheta;
    const G4double cos2Theta = cosTheta * cosTheta;

    ++event.num_encounter_total;
    event.sum_cos_theta_encounter += cosTheta;
    event.sum_one_minus_cos_theta_encounter += oneMinusCosTheta;
    event.sum_cos2_theta_encounter += cos2Theta;
    ++event.phase_function_histogram[PhaseFunctionBin(cosTheta)];

    ++event.num_particle_scatter;
    event.sum_cos_theta_particle += cosTheta;

    if (prePhase == DetectorConstruction::Phase::BN)
    {
      ++event.num_encounter_BN;
      event.sum_cos_theta_encounter_BN += cosTheta;
      event.sum_one_minus_cos_theta_encounter_BN += oneMinusCosTheta;
      event.sum_cos2_theta_encounter_BN += cos2Theta;

      ++event.num_particle_scatter_BN;
      event.sum_cos_theta_particle_BN += cosTheta;
    }
    else if (prePhase == DetectorConstruction::Phase::ZnS)
    {
      ++event.num_encounter_ZnS;
      event.sum_cos_theta_encounter_ZnS += cosTheta;
      event.sum_one_minus_cos_theta_encounter_ZnS += oneMinusCosTheta;
      event.sum_cos2_theta_encounter_ZnS += cos2Theta;

      ++event.num_particle_scatter_ZnS;
      event.sum_cos_theta_particle_ZnS += cosTheta;
    }

    event.in_particle_segment = false;
    event.particle_segment_phase.clear();
  }

  if (!isOuterRveExit && processName != "OpAbsorption")
  {
    const G4double thetaDeg =
        AngleDeg(preDir, postDir);
    if (thetaDeg > fConfig->stageD_theta_threshold_deg)
    {
      ++event.num_real_scatter;
      const G4double cosTheta = std::clamp(
          preDir.dot(postDir),
          -1.0, 1.0);
      event.sum_cos_theta += cosTheta;
      if (isMaterialBoundary)
      {
        ++event.num_boundary_scatter;
        event.sum_cos_theta_boundary += cosTheta;
        if (prePhase == DetectorConstruction::Phase::BN)
        {
          ++event.num_boundary_scatter_BN;
          event.sum_cos_theta_boundary_BN += cosTheta;
        }
        else if (prePhase == DetectorConstruction::Phase::ZnS)
        {
          ++event.num_boundary_scatter_ZnS;
          event.sum_cos_theta_boundary_ZnS += cosTheta;
        }
      }
      else
      {
        ++event.num_bulk_scatter;
        event.sum_cos_theta_bulk += cosTheta;
      }
    }
  }

  if (detector != nullptr && HandleBoundaryReentry(step, track, detector))
    return;

  HandleLimitKills(step, track);
}
