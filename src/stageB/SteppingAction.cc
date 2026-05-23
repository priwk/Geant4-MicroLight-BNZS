#include "SteppingAction.hh"

#include "EventAction.hh"
#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"

#include <cmath>
#include <fstream>
#include <string>

// --------------------------------------------------------------------
// helpers
namespace
{
    G4bool StartsWith(const G4String &s, const char *prefix)
    {
        const std::string value = s;
        const std::string p = prefix;
        return value.rfind(p, 0) == 0;
    }

    G4bool IsTrackedHeavyParticle(const G4Track *track)
    {
        const auto *def = track->GetDefinition();
        if (!def)
            return false;

        if (def->GetParticleName() == "alpha")
        {
            return true;
        }

        // Li7 ion
        if (def->GetParticleType() == "nucleus" &&
            def->GetAtomicNumber() == 3 &&
            def->GetAtomicMass() == 7)
        {
            return true;
        }

        return false;
    }

    std::string ParticleLabel(const G4Track *track)
    {
        const auto *def = track->GetDefinition();
        if (!def)
            return "unknown";

        if (def->GetParticleName() == "alpha")
        {
            return "alpha";
        }

        if (def->GetParticleType() == "nucleus" &&
            def->GetAtomicNumber() == 3 &&
            def->GetAtomicMass() == 7)
        {
            return "Li7";
        }

        return def->GetParticleName();
    }

    std::string PhaseLabel(const G4VPhysicalVolume *pv)
    {
        if (!pv)
            return "outside";

        const auto *lv = pv->GetLogicalVolume();
        if (!lv)
            return "outside";

        const auto &lvName = lv->GetName();

        if (lvName == "BN_LV" || StartsWith(lvName, "BN_ClipLV_"))
            return "BN";
        if (lvName == "ZnS_LV" || StartsWith(lvName, "ZnS_ClipLV_"))
            return "ZnS";
        if (lvName == "MatrixLV")
            return "binder_void";
        if (lvName == "WorldLV")
            return "outside";

        return "other";
    }

    G4bool IsOutsidePhase(const std::string &phase)
    {
        return (phase == "outside");
    }

    G4double PatchHalfXYUm(const DetectorConstruction *det)
    {
        return 0.5 * det->GetPatchXYUm();
    }

    G4double PatchHalfZUm(const DetectorConstruction *det)
    {
        return 0.5 * det->GetEffectiveLocalThickness() / um;
    }

    std::string ExitFaceLabel(const G4ThreeVector &position, const DetectorConstruction *det)
    {
        const G4double tolUm = 1.0e-3;
        const G4double xUm = position.x() / um;
        const G4double yUm = position.y() / um;
        const G4double zUm = position.z() / um;
        const G4double halfXYUm = PatchHalfXYUm(det);
        const G4double halfZUm = PatchHalfZUm(det);

        if (zUm >= halfZUm - tolUm)
            return "+z";
        if (zUm <= -halfZUm + tolUm)
            return "-z";
        if (xUm >= halfXYUm - tolUm)
            return "+x";
        if (xUm <= -halfXYUm + tolUm)
            return "-x";
        if (yUm >= halfXYUm - tolUm)
            return "+y";
        if (yUm <= -halfXYUm + tolUm)
            return "-y";
        return "unknown";
    }

    RunAction::BoundaryExitClass ClassifyBoundaryExit(
        const std::string &surfaceMode,
        const std::string &exitFace)
    {
        if (surfaceMode == "front_surface" && exitFace == "+z")
        {
            return RunAction::BoundaryExitClass::PhysicalSurfaceExit;
        }
        if (surfaceMode == "back_surface" && exitFace == "-z")
        {
            return RunAction::BoundaryExitClass::PhysicalSurfaceExit;
        }
        return RunAction::BoundaryExitClass::UnexpectedArtificialExit;
    }
}

// --------------------------------------------------------------------

SteppingAction::SteppingAction(EventAction *eventAction,
                               const PrimaryGeneratorAction *primaryAction)
    : G4UserSteppingAction(),
      fEventAction(eventAction),
      fPrimaryAction(primaryAction)
{
}

// --------------------------------------------------------------------

SteppingAction::~SteppingAction() = default;

// --------------------------------------------------------------------

void SteppingAction::UserSteppingAction(const G4Step *step)
{
    if (!step)
        return;

    auto *track = step->GetTrack();
    if (!track)
        return;

    if (!IsTrackedHeavyParticle(track))
        return;

    const auto *prePoint = step->GetPreStepPoint();
    const auto *postPoint = step->GetPostStepPoint();
    if (!prePoint || !postPoint)
        return;

    const auto *prePV = prePoint->GetPhysicalVolume();
    const auto *postPV = postPoint->GetPhysicalVolume();

    const std::string phasePre = PhaseLabel(prePV);
    const std::string phasePost = PhaseLabel(postPV);

    // Do not keep tracking into world vacuum after leaving the patch.
    // But record the boundary-crossing step itself.
    const G4ThreeVector &xPre = prePoint->GetPosition();
    const G4ThreeVector &xPost = postPoint->GetPosition();

    const G4double stepLen = step->GetStepLength();
    const G4double edep = step->GetTotalEnergyDeposit();
    const G4double ekinPre = prePoint->GetKineticEnergy();
    const G4double ekinPost = postPoint->GetKineticEnergy();

    if (fEventAction)
    {
        fEventAction->AddEdep(edep);
    }

    auto *runAction = fEventAction ? fEventAction->GetRunAction() : nullptr;
    if (runAction && fPrimaryAction)
    {
        runAction->SwitchOutputCsvForInputPath(fPrimaryAction->GetLoadedInputFile());
    }

    if (runAction && fPrimaryAction && runAction->IsFullMode() && runAction->IsFullStepCsvOpen())
    {
        std::ofstream &csv = runAction->GetFullStepCsv();

        const auto &rec = fPrimaryAction->GetCurrentRecord();
        const auto &capturePos = fPrimaryAction->GetCurrentLocalCapturePosition();
        const auto &bnCenter = fPrimaryAction->GetCurrentSelectedBNCenter();
        std::string placementFile = "unknown";

        const auto *det = dynamic_cast<const DetectorConstruction *>(
            G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        if (det)
        {
            placementFile = det->GetLoadedPlacementFileForRecord();
        }

        csv
            << fPrimaryAction->MakeCurrentPhysicalEventUid() << ","
            << rec.eventID << ","
            << rec.thickness_um << ","
            << rec.bn_wt << ","
            << rec.zns_wt << ","
            << rec.capture_x_um << ","
            << rec.capture_y_um << ","
            << rec.corr_x_um << ","
            << rec.corr_y_um << ","
            << rec.depth_um << ","
            << placementFile << ","
            << capturePos.x() / um << ","
            << capturePos.y() / um << ","
            << capturePos.z() / um << ","
            << fPrimaryAction->GetCurrentSurfaceMode() << ","
            << fPrimaryAction->GetCurrentTargetLocalZ() / um << ","
            << fPrimaryAction->GetCurrentUsedLocalZ() / um << ","
            << bnCenter.x() / um << ","
            << bnCenter.y() / um << ","
            << bnCenter.z() / um << ","
            << fPrimaryAction->GetCurrentAlphaLiReplayIndex() << ","
            << fPrimaryAction->GetCurrentAlphaLiReplayCount() << ","
            << track->GetTrackID() << ","
            << track->GetCurrentStepNumber() << ","
            << ParticleLabel(track) << ","
            << phasePre << ","
            << phasePost << ","
            << xPre.x() / um << ","
            << xPre.y() / um << ","
            << xPre.z() / um << ","
            << xPost.x() / um << ","
            << xPost.y() / um << ","
            << xPost.z() / um << ","
            << stepLen / um << ","
            << edep / keV << ","
            << ekinPre / keV << ","
            << ekinPost / keV << ","
            << fPrimaryAction->MakeCurrentSourceEventUid() << ","
            << rec.record_index << ","
            << fPrimaryAction->GetCurrentTrajectoryWeight()
            << "\n";
    }
    else if (runAction && fEventAction && runAction->IsSlimMode() &&
             runAction->IsSlimTrackCsvOpen() &&
             phasePre == "ZnS" && stepLen > 0.0)
    {
        std::ofstream &csv = runAction->GetSlimTrackCsv();
        const auto anchor = fEventAction->MakeCurrentCaptureAnchorRow();

        csv
            << anchor.physical_event_uid << ","
            << anchor.source_event_uid << ","
            << anchor.eventID << ","
            << anchor.record_index << ","
            << track->GetTrackID() << ","
            << track->GetCurrentStepNumber() << ","
            << ParticleLabel(track) << ","
            << phasePost << ","
            << xPre.x() / um << ","
            << xPre.y() / um << ","
            << xPre.z() / um << ","
            << xPost.x() / um << ","
            << xPost.y() / um << ","
            << xPost.z() / um << ","
            << stepLen / um << ","
            << edep / keV << ","
            << ekinPre / keV << ","
            << ekinPost / keV << ","
            << anchor.alphali_replay_index << ","
            << anchor.alphali_replay_count << ","
            << anchor.trajectory_weight
            << "\n";
    }

    // Kill track once it exits the microstructure into world/outside
    if (!IsOutsidePhase(phasePre) && IsOutsidePhase(phasePost))
    {
        if (runAction && fEventAction)
        {
            const auto *det = dynamic_cast<const DetectorConstruction *>(
                G4RunManager::GetRunManager()->GetUserDetectorConstruction());
            if (det)
            {
                const auto anchor = fEventAction->MakeCurrentCaptureAnchorRow();
                const std::string exitFace = ExitFaceLabel(xPost, det);
                runAction->RecordBoundaryExit(
                    anchor,
                    ParticleLabel(track),
                    ClassifyBoundaryExit(anchor.surface_mode, exitFace),
                    ekinPost / keV);
            }
        }
        track->SetTrackStatus(fStopAndKill);
    }
}
