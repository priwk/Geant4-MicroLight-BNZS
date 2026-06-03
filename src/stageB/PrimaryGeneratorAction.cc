#include "PrimaryGeneratorAction.hh"

#include "AnalysisConfig.hh"
#include "DetectorConstruction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Exception.hh"
#include "G4ExceptionSeverity.hh"

#include "Randomize.hh"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdlib>
// #include <dirent.h>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <numeric>

// --------------------------------------------------------------------
// small helpers in anonymous namespace
namespace
{
  std::vector<std::string> SplitFlexible(const std::string &line)
  {
    std::string s = line;
    for (char &c : s)
    {
      if (c == ',')
        c = '\t';
    }

    std::vector<std::string> tokens;
    std::stringstream ss(s);
    std::string token;
    while (std::getline(ss, token, '\t'))
    {
      if (!token.empty() && token.back() == '\r')
        token.pop_back();
      tokens.push_back(token);
    }
    return tokens;
  }

  std::string Trim(const std::string &value)
  {
    const auto begin = value.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos)
      return "";
    const auto end = value.find_last_not_of(" \t\r\n");
    return value.substr(begin, end - begin + 1);
  }

  G4int PositiveIntFromEnv(const char *name, G4int fallback)
  {
    const char *raw = std::getenv(name);
    if (raw == nullptr)
      return fallback;

    try
    {
      const G4int value = std::stoi(Trim(raw));
      return (value > 0) ? value : fallback;
    }
    catch (...)
    {
      return fallback;
    }
  }

  G4ThreeVector RandomUnitVector()
  {
    const G4double cosTheta = 2.0 * G4UniformRand() - 1.0;
    const G4double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
    const G4double phi = twopi * G4UniformRand();

    return G4ThreeVector(
        sinTheta * std::cos(phi),
        sinTheta * std::sin(phi),
        cosTheta);
  }

  G4bool EndsWith(const std::string &s, const std::string &suffix)
  {
    if (s.size() < suffix.size())
      return false;
    return s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
  }

  std::filesystem::path MakeInputRootDirectoryPath()
  {
    return std::filesystem::current_path().parent_path() / "Input";
  }

  std::filesystem::path MakeCaptureInputRootDirectoryPath()
  {
    return MakeInputRootDirectoryPath() / "neutron_capture_positions";
  }

  std::filesystem::path MakeStageACaptureInputRootDirectoryPath()
  {
    return MakeInputRootDirectoryPath() / "stageA";
  }

  G4bool TryParseRatioFolderName(const std::string &name,
                                 G4double &bnWt,
                                 G4double &znsWt)
  {
    const std::size_t dashPos = name.find('-');
    if (dashPos == std::string::npos)
      return false;

    const std::string lhs = name.substr(0, dashPos);
    const std::string rhs = name.substr(dashPos + 1);

    if (lhs.empty() || rhs.empty())
      return false;

    try
    {
      bnWt = std::stod(lhs);
      znsWt = std::stod(rhs);
    }
    catch (...)
    {
      return false;
    }

    return (bnWt > 0.0 && znsWt > 0.0);
  }

  std::vector<std::string> CollectMatchingCsvFilesInDirectory(const std::filesystem::path &dir)
  {
    namespace fs = std::filesystem;
    std::vector<std::string> matches;

    if (!fs::exists(dir) || !fs::is_directory(dir))
      return matches;

    for (const auto &entry : fs::directory_iterator(dir))
    {
      if (!entry.is_regular_file())
        continue;

      const std::string name = entry.path().filename().string();
      if (EndsWith(name, "_neutron_capture_positions.csv"))
      {
        matches.push_back(entry.path().string());
      }
    }

    auto thicknessTag = [](const std::string &path) -> G4double
    {
      const auto fileName = std::filesystem::path(path).filename().string();
      const std::string key = "_neutron_capture_positions.csv";
      const auto keyPos = fileName.find(key);
      if (keyPos == std::string::npos)
        return DBL_MAX;

      try
      {
        return std::stod(fileName.substr(0, keyPos));
      }
      catch (...)
      {
        return DBL_MAX;
      }
    };

    std::sort(matches.begin(), matches.end(),
              [&](const std::string &a, const std::string &b)
              {
                const G4double ta = thicknessTag(a);
                const G4double tb = thicknessTag(b);
                if (std::abs(ta - tb) > 1.0e-12)
                  return ta < tb;
                return a < b;
              });
    return matches;
  }

  std::string WeightPartToTagString(G4double v)
  {
    const G4double rounded = std::round(v);
    if (std::fabs(v - rounded) < 1.0e-9)
    {
      std::ostringstream oss;
      oss << static_cast<long long>(rounded);
      return oss.str();
    }

    std::ostringstream oss;
    oss << v;
    return oss.str();
  }

  std::string MakeRatioFolderName(G4double bnWt, G4double znsWt)
  {
    return WeightPartToTagString(bnWt) + "-" + WeightPartToTagString(znsWt);
  }

  const G4double kXYMargin = 7.0 * um;
  const G4double kZBulkMargin = 7.45 * um;

  G4bool IsInsideSafeXY(const G4ThreeVector &p, const DetectorConstruction *det)
  {
    const G4double halfXY = 0.5 * det->GetPatchXYUm() * um;
    return (std::abs(p.x()) <= halfXY - kXYMargin &&
            std::abs(p.y()) <= halfXY - kXYMargin);
  }

  G4bool IsInsideBulkZWindow(const G4ThreeVector &p, const DetectorConstruction *det)
  {
    const G4double localT = det->GetEffectiveLocalThickness();
    return (p.z() >= -0.5 * localT + kZBulkMargin &&
            p.z() <= +0.5 * localT - kZBulkMargin);
  }

  G4bool IsBulkSphereFullyInSafeWindow(const G4ThreeVector &center,
                                       G4double radius,
                                       const DetectorConstruction *det)
  {
    const G4double halfXY = 0.5 * det->GetPatchXYUm() * um;
    const G4double localT = det->GetEffectiveLocalThickness();

    return (std::abs(center.x()) + radius <= halfXY - kXYMargin &&
            std::abs(center.y()) + radius <= halfXY - kXYMargin &&
            center.z() - radius >= -0.5 * localT + kZBulkMargin &&
            center.z() + radius <= +0.5 * localT - kZBulkMargin);
  }

  G4bool IsBulkPointSafe(const G4ThreeVector &p, const DetectorConstruction *det)
  {
    return IsInsideSafeXY(p, det) && IsInsideBulkZWindow(p, det);
  }

  G4bool DoesSphereIntersectBulkSafeWindow(const G4ThreeVector &center,
                                           G4double radius,
                                           const DetectorConstruction *det)
  {
    const G4double safeHalfXY = 0.5 * det->GetPatchXYUm() * um - kXYMargin;
    const G4double localT = det->GetEffectiveLocalThickness();
    const G4double zMin = -0.5 * localT + kZBulkMargin;
    const G4double zMax = +0.5 * localT - kZBulkMargin;

    const G4double closestX = std::min(std::max(center.x(), -safeHalfXY), safeHalfXY);
    const G4double closestY = std::min(std::max(center.y(), -safeHalfXY), safeHalfXY);
    const G4double closestZ = std::min(std::max(center.z(), zMin), zMax);

    const G4double dx = center.x() - closestX;
    const G4double dy = center.y() - closestY;
    const G4double dz = center.z() - closestZ;

    return (dx * dx + dy * dy + dz * dz <= radius * radius);
  }

  G4bool IsSurfaceSliceSafelyInsideXY(const G4ThreeVector &center,
                                      G4double diskRadius,
                                      const DetectorConstruction *det)
  {
    const G4double halfXY = 0.5 * det->GetPatchXYUm() * um;
    return (std::abs(center.x()) + diskRadius <= halfXY - kXYMargin &&
            std::abs(center.y()) + diskRadius <= halfXY - kXYMargin);
  }

  G4bool DoesSurfaceSliceIntersectSafeXY(const G4ThreeVector &center,
                                         G4double diskRadius,
                                         const DetectorConstruction *det)
  {
    const G4double safeHalfXY = 0.5 * det->GetPatchXYUm() * um - kXYMargin;
    const G4double dx = std::max(0.0, std::abs(center.x()) - safeHalfXY);
    const G4double dy = std::max(0.0, std::abs(center.y()) - safeHalfXY);
    return (dx * dx + dy * dy <= diskRadius * diskRadius);
  }
}

// --------------------------------------------------------------------

PrimaryGeneratorAction::PrimaryGeneratorAction(AnalysisConfig *config)
    : G4VUserPrimaryGeneratorAction(),
      fConfig(config),
      fParticleGun(nullptr),
      fInputFiles(),
      fCurrentFileIndex(0),
      fCurrentInputStream(),
      fCurrentInputFile(""),
      fCurrentRecordInputFile(""),
      fCurrentHeaderIndex(),
      fCurrentInputRecordCounter(0),
      fFirstRecordForGeometry(),
      fHasFirstRecordForGeometry(false),
      fNoMoreInput(false),
      fTotalStreamedRecords(0),
      fAlphaLiReplayPerCapture(ReadAlphaLiReplayPerCapture()),
      fCurrentAlphaLiReplayIndex(0),
      fRemainingReplaysForCurrentCapture(0),
      fCurrentRecord(),
      fCurrentLocalCapturePosition(),
      fCurrentSelectedBNCenter(),
      fCurrentSurfaceMode(""),
      fCurrentTargetLocalZ(0.0),
      fCurrentUsedLocalZ(0.0)
{
  fParticleGun = new G4ParticleGun(1);

  InitializeInputStreaming();
  ConfigureDetectorFromInput();

  // temporary default, will be overwritten event-by-event
  auto *particleTable = G4ParticleTable::GetParticleTable();
  auto *alpha = particleTable->FindParticle("alpha");
  fParticleGun->SetParticleDefinition(alpha);
  fParticleGun->SetParticleEnergy(1.0 * MeV);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
  fParticleGun->SetParticlePosition(G4ThreeVector());
}

// --------------------------------------------------------------------

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

// --------------------------------------------------------------------

std::vector<std::string> PrimaryGeneratorAction::FindInputCsvFiles() const
{
  // Priority 0: AnalysisConfig -> single explicit file
  if (fConfig != nullptr && !fConfig->captureCsvPath.empty())
  {
    return {fConfig->captureCsvPath};
  }

  // Priority 1: environment variable -> single explicit file
  const char *envPath = std::getenv("BNZS_INPUT_CSV");
  if (envPath && std::string(envPath).size() > 0)
  {
    return {std::string(envPath)};
  }

  namespace fs = std::filesystem;

  // Priority 2: AnalysisConfig/env -> explicit ratio directory
  std::string explicitDir;
  if (fConfig != nullptr && !fConfig->captureInputDir.empty())
  {
    explicitDir = fConfig->captureInputDir;
  }
  else
  {
    const char *envDir = std::getenv("BNZS_INPUT_DIR");
    if (envDir && std::string(envDir).size() > 0)
    {
      explicitDir = envDir;
    }
  }

  if (!explicitDir.empty())
  {
    auto matches = CollectMatchingCsvFilesInDirectory(explicitDir);
    if (matches.empty())
    {
      G4Exception("PrimaryGeneratorAction::FindInputCsvFiles",
                  "BNZS001", FatalException,
                  ("No *_neutron_capture_positions.csv found in explicit Stage B input directory: " + explicitDir).c_str());
      return {};
    }
    return matches;
  }

  // Priority 3: preferred new layout, selected by current BN:ZnS ratio
  // ../Input/stageA/<ratio>/neutron_capture_positions/*.csv
  // Backward-compatible fallback:
  // ../Input/neutron_capture_positions/<ratio>/*.csv
  const std::vector<fs::path> captureRoots = {
      MakeStageACaptureInputRootDirectoryPath(),
      MakeCaptureInputRootDirectoryPath(),
  };

  for (const auto &captureRoot : captureRoots)
  {
    if (!fs::exists(captureRoot) || !fs::is_directory(captureRoot))
    {
      continue;
    }

    if (fConfig != nullptr)
    {
      fs::path ratioDir = captureRoot / MakeRatioFolderName(fConfig->bnWt, fConfig->znsWt);
      if (captureRoot.filename() == "stageA")
      {
        ratioDir /= "neutron_capture_positions";
      }
      auto matches = CollectMatchingCsvFilesInDirectory(ratioDir);
      if (!matches.empty())
      {
        return matches;
      }
    }

    std::vector<fs::path> candidateRatioDirs;

    for (const auto &entry : fs::directory_iterator(captureRoot))
    {
      if (!entry.is_directory())
        continue;

      G4double bnWt = 0.0, znsWt = 0.0;
      const std::string dirName = entry.path().filename().string();

      if (!TryParseRatioFolderName(dirName, bnWt, znsWt))
        continue;

      fs::path ratioDir = entry.path();
      if (captureRoot.filename() == "stageA")
      {
        ratioDir /= "neutron_capture_positions";
      }

      auto matches = CollectMatchingCsvFilesInDirectory(ratioDir);
      if (!matches.empty())
      {
        candidateRatioDirs.push_back(ratioDir);
      }
    }

    if (candidateRatioDirs.size() == 1)
    {
      return CollectMatchingCsvFilesInDirectory(candidateRatioDirs.front());
    }

    if (candidateRatioDirs.size() > 1)
    {
      std::ostringstream oss;
      oss << "Multiple non-empty ratio folders found under "
          << captureRoot.string()
          << ". Set /cfg/setWeightRatio, /cfg/setCaptureDir, BNZS_INPUT_DIR, or BNZS_INPUT_CSV for one run.";
      G4Exception("PrimaryGeneratorAction::FindInputCsvFiles",
                  "BNZS001", FatalException,
                  oss.str().c_str());
      return {};
    }
  }

  // Priority 4: backward compatibility
  // old flat layout: ../Input/*_neutron_capture_positions.csv
  const fs::path inputDir = MakeInputRootDirectoryPath();
  auto matches = CollectMatchingCsvFilesInDirectory(inputDir);

  if (matches.empty())
  {
    G4Exception("PrimaryGeneratorAction::FindInputCsvFiles",
                "BNZS002", FatalException,
                ("No *_neutron_capture_positions.csv found in " +
                 MakeStageACaptureInputRootDirectoryPath().string() + ", " +
                 MakeCaptureInputRootDirectoryPath().string() + ", or " +
                 inputDir.string())
                    .c_str());
    return {};
  }

  return matches;
}

G4bool PrimaryGeneratorAction::ReadHeaderIndex(std::istream &input, HeaderIndex &headerIndex) const
{
  std::string headerLine;
  if (!std::getline(input, headerLine))
    return false;

  headerIndex.clear();
  const auto headers = SplitFlexible(headerLine);
  for (std::size_t i = 0; i < headers.size(); ++i)
  {
    headerIndex[Trim(headers[i])] = i;
  }

  return !headerIndex.empty();
}

G4bool PrimaryGeneratorAction::ParseOneRecordLine(
    const std::string &line,
    const HeaderIndex &headerIndex,
    G4int fallbackRecordIndex,
    CaptureRecord &rec) const
{
  const auto tokens = SplitFlexible(line);

  auto field = [&](const std::string &name) -> std::string
  {
    const auto it = headerIndex.find(name);
    if (it == headerIndex.end() || it->second >= tokens.size())
      return "";
    return Trim(tokens[it->second]);
  };

  const std::string eventID = field("eventID");
  const std::string thickness = field("thickness_um");
  const std::string bnWt = field("bn_wt");
  const std::string znsWt = field("zns_wt");
  const std::string captureX = field("capture_x_um");
  const std::string captureY = field("capture_y_um");
  std::string sourceX = field("source_x_um");
  std::string sourceY = field("source_y_um");
  if (sourceX.empty())
    sourceX = field("corr_x_um");
  if (sourceY.empty())
    sourceY = field("corr_y_um");
  const std::string depth = field("depth_um");

  if (eventID.empty() || thickness.empty() || bnWt.empty() || znsWt.empty() ||
      captureX.empty() || captureY.empty() || sourceX.empty() || sourceY.empty() ||
      depth.empty())
  {
    return false;
  }

  try
  {
    rec.eventID = std::stoi(eventID);
    rec.thickness_um = std::stod(thickness);
    rec.bn_wt = std::stod(bnWt);
    rec.zns_wt = std::stod(znsWt);
    rec.capture_x_um = std::stod(captureX);
    rec.capture_y_um = std::stod(captureY);
    rec.source_x_um = std::stod(sourceX);
    rec.source_y_um = std::stod(sourceY);
    rec.depth_um = std::stod(depth);
  }
  catch (...)
  {
    return false;
  }

  const std::string recordIndex = field("record_index");
  if (!recordIndex.empty())
  {
    try
    {
      rec.record_index = std::stoi(recordIndex);
    }
    catch (...)
    {
      rec.record_index = fallbackRecordIndex;
    }
  }
  else
  {
    rec.record_index = fallbackRecordIndex;
  }

  return true;
}

// --------------------------------------------------------------------

G4bool PrimaryGeneratorAction::ReadFirstValidRecordFromFile(
    const std::string &path, CaptureRecord &rec) const
{
  std::ifstream fin(path.c_str());
  if (!fin)
    return false;

  HeaderIndex headerIndex;
  if (!ReadHeaderIndex(fin, headerIndex))
    return false;

  std::string line;
  G4int fallbackRecordIndex = 0;
  while (std::getline(fin, line))
  {
    if (line.empty())
      continue;

    if (ParseOneRecordLine(line, headerIndex, fallbackRecordIndex, rec))
    {
      return true;
    }

    ++fallbackRecordIndex;
  }

  return false;
}

G4bool PrimaryGeneratorAction::OpenNextInputFile()
{
  if (fCurrentInputStream.is_open())
  {
    fCurrentInputStream.close();
  }

  while (fCurrentFileIndex < fInputFiles.size())
  {
    fCurrentInputFile = fInputFiles[fCurrentFileIndex++];
    fCurrentInputStream.open(fCurrentInputFile.c_str());

    if (!fCurrentInputStream)
    {
      continue;
    }

    if (!ReadHeaderIndex(fCurrentInputStream, fCurrentHeaderIndex))
    {
      fCurrentInputStream.close();
      continue;
    }

    fCurrentInputRecordCounter = 0;

    G4cout << "[PrimaryGeneratorAction] Open input CSV: "
           << fCurrentInputFile << G4endl;
    return true;
  }

  fNoMoreInput = true;
  return false;
}

G4bool PrimaryGeneratorAction::ReadNextRecord(CaptureRecord &rec)
{
  std::string line;

  while (true)
  {
    if (!fCurrentInputStream.is_open())
    {
      if (!OpenNextInputFile())
      {
        return false;
      }
    }

    while (std::getline(fCurrentInputStream, line))
    {
      if (line.empty())
        continue;

      if (ParseOneRecordLine(
              line,
              fCurrentHeaderIndex,
              fCurrentInputRecordCounter,
              rec))
      {
        fCurrentRecordInputFile = fCurrentInputFile;
        ++fTotalStreamedRecords;
        ++fCurrentInputRecordCounter;
        return true;
      }

      ++fCurrentInputRecordCounter;
    }

    fCurrentInputStream.close();
  }
}

void PrimaryGeneratorAction::InitializeInputStreaming()
{
  const auto candidateInputFiles = FindInputCsvFiles();
  fInputFiles.clear();
  fCurrentFileIndex = 0;
  fNoMoreInput = false;
  fTotalStreamedRecords = 0;
  fHasFirstRecordForGeometry = false;
  fCurrentInputFile.clear();
  fCurrentRecordInputFile.clear();
  fCurrentHeaderIndex.clear();
  fCurrentInputRecordCounter = 0;
  fCurrentAlphaLiReplayIndex = 0;
  fRemainingReplaysForCurrentCapture = 0;

  if (candidateInputFiles.empty())
  {
    G4Exception("PrimaryGeneratorAction::InitializeInputStreaming",
                "BNZS003", FatalException,
                "No input CSV files found.");
    return;
  }

  const auto *det = dynamic_cast<const DetectorConstruction *>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  const G4double localT_um = det ? det->GetEffectiveLocalThickness() / um : 30.0;

  for (const auto &path : candidateInputFiles)
  {
    CaptureRecord rec;
    if (!ReadFirstValidRecordFromFile(path, rec))
    {
      G4Exception("PrimaryGeneratorAction::InitializeInputStreaming",
                  "BNZS005", FatalException,
                  ("Failed to read first valid record from: " + path).c_str());
      return;
    }

    if (!IsInputThicknessCompatible(rec.thickness_um, localT_um))
    {
      G4cout << "[PrimaryGeneratorAction] Skip input CSV thinner than local patch: "
             << path
             << " (input thickness = " << rec.thickness_um
             << " um, local patch thickness = " << localT_um << " um)"
             << G4endl;
      continue;
    }

    fInputFiles.push_back(path);
  }

  if (fInputFiles.empty())
  {
    G4Exception("PrimaryGeneratorAction::InitializeInputStreaming",
                "BNZS003", FatalException,
                "No Stage B input CSV files are compatible with the local microstructure thickness.");
    return;
  }

  CaptureRecord refRec;
  if (!ReadFirstValidRecordFromFile(fInputFiles.front(), refRec))
  {
    G4Exception("PrimaryGeneratorAction::InitializeInputStreaming",
                "BNZS004", FatalException,
                ("Failed to read first valid record from: " + fInputFiles.front()).c_str());
    return;
  }

  fFirstRecordForGeometry = refRec;
  fHasFirstRecordForGeometry = true;

  for (const auto &path : fInputFiles)
  {
    CaptureRecord rec;
    if (!ReadFirstValidRecordFromFile(path, rec))
    {
      G4Exception("PrimaryGeneratorAction::InitializeInputStreaming",
                  "BNZS005", FatalException,
                  ("Failed to read first valid record from: " + path).c_str());
      return;
    }

    if (std::abs(rec.bn_wt - refRec.bn_wt) > 1e-12 ||
        std::abs(rec.zns_wt - refRec.zns_wt) > 1e-12)
    {
      G4Exception("PrimaryGeneratorAction::InitializeInputStreaming",
                  "BNZS006", FatalException,
                  ("Mixed BN:ZnS ratios are not allowed across input files: " + path).c_str());
      return;
    }
  }

  if (!OpenNextInputFile())
  {
    G4Exception("PrimaryGeneratorAction::InitializeInputStreaming",
                "BNZS008", FatalException,
                "Failed to open the first input CSV.");
    return;
  }

  G4cout
      << "\n[PrimaryGeneratorAction] Streaming input enabled"
      << "\n  files = " << fInputFiles.size()
      << "\n  fixed BN:ZnS wt = " << refRec.bn_wt << " : " << refRec.zns_wt
      << "\n  local geometry thickness = " << localT_um << " um"
      << G4endl;
}

// --------------------------------------------------------------------

void PrimaryGeneratorAction::ConfigureDetectorFromInput()
{
  const auto *detConst = dynamic_cast<const DetectorConstruction *>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  auto *det = const_cast<DetectorConstruction *>(detConst);

  if (!det || !fHasFirstRecordForGeometry)
  {
    G4Exception("PrimaryGeneratorAction::ConfigureDetectorFromInput",
                "BNZS009", FatalException,
                "DetectorConstruction is not available or first input record missing.");
    return;
  }

  // Keep microstructure geometry fixed.
  // Only set composition once from the first input record.
  det->SetWeightRatio(fFirstRecordForGeometry.bn_wt,
                      fFirstRecordForGeometry.zns_wt);
}

// --------------------------------------------------------------------

G4bool PrimaryGeneratorAction::IsInputThicknessCompatible(
    G4double thickness_um, G4double localT_um) const
{
  const G4double tol = 1.0e-12;

  // 默认新行为：允许 thickness == local patch thickness
  const G4bool allowEqual =
      (fConfig == nullptr) ? true : fConfig->allowThicknessEqualLocalPatch;

  if (allowEqual)
  {
    return (thickness_um + tol >= localT_um);
  }

  return (thickness_um > localT_um + tol);
}

// --------------------------------------------------------------------

G4int PrimaryGeneratorAction::ReadAlphaLiReplayPerCapture() const
{
  return PositiveIntFromEnv("BNZS_ALPHALI_REPLAY_PER_CAPTURE", 1);
}

G4bool PrimaryGeneratorAction::PrepareCurrentCaptureReplayState()
{
  CaptureRecord rec;
  if (!ReadNextRecord(rec))
  {
    return false;
  }

  fCurrentRecord = rec;
  fCurrentSurfaceMode = DetermineSurfaceMode(fCurrentRecord);
  fCurrentTargetLocalZ = DetermineTargetLocalZ(fCurrentRecord, fCurrentSurfaceMode);

  const auto *det = dynamic_cast<const DetectorConstruction *>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  G4ThreeVector chosenCenter;
  G4double usedZ = 0.0;
  G4bool usedFallback = false;

  if (fCurrentSurfaceMode == "bulk")
  {
    if (!SampleBulkCapturePoint(chosenCenter, fCurrentLocalCapturePosition, usedFallback))
    {
      G4Exception("PrimaryGeneratorAction::PrepareCurrentCaptureReplayState",
                  "BNZS008", FatalException,
                  "Failed to sample a safe bulk BN capture point.");
      return false;
    }

    usedZ = fCurrentLocalCapturePosition.z();
    fCurrentTargetLocalZ = usedZ;
  }
  else
  {
    if (!SelectBNSphereForTargetZ(
            fCurrentTargetLocalZ,
            chosenCenter,
            usedZ,
            usedFallback))
    {
      G4Exception("PrimaryGeneratorAction::PrepareCurrentCaptureReplayState",
                  "BNZS008", FatalException,
                  "Failed to find any BN sphere for current surface capture event.");
      return false;
    }

    if (!SampleSafePointInSphereSlice(
            chosenCenter,
            usedZ,
            det->GetBNRadius(),
            fCurrentLocalCapturePosition))
    {
      G4Exception("PrimaryGeneratorAction::PrepareCurrentCaptureReplayState",
                  "BNZS010", FatalException,
                  "Failed to sample a surface capture point inside the XY safety window.");
      return false;
    }
  }

  fCurrentSelectedBNCenter = chosenCenter;
  fCurrentUsedLocalZ = usedZ;
  fCurrentAlphaLiReplayIndex = 0;
  fRemainingReplaysForCurrentCapture = fAlphaLiReplayPerCapture;

  const G4int geantEventID = G4RunManager::GetRunManager()->GetCurrentEvent()
                                 ? G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()
                                 : -1;
  if (geantEventID >= 0 && geantEventID < 5)
  {
    G4cout
        << "\n[PrimaryGeneratorAction] Event " << geantEventID
        << "\n  input eventID      = " << fCurrentRecord.eventID
        << "\n  record index       = " << fCurrentRecord.record_index
        << "\n  mode               = " << fCurrentSurfaceMode
        << "\n  target z           = " << fCurrentTargetLocalZ / um << " um"
        << "\n  used z             = " << fCurrentUsedLocalZ / um << " um"
        << "\n  selected BN center = (" << fCurrentSelectedBNCenter.x() / um
        << ", " << fCurrentSelectedBNCenter.y() / um
        << ", " << fCurrentSelectedBNCenter.z() / um << ") um"
        << "\n  capture point      = (" << fCurrentLocalCapturePosition.x() / um
        << ", " << fCurrentLocalCapturePosition.y() / um
        << ", " << fCurrentLocalCapturePosition.z() / um << ") um"
        << "\n  replay count       = " << fAlphaLiReplayPerCapture
        << "\n  fallback used      = " << (usedFallback ? "yes" : "no")
        << G4endl;
  }

  return true;
}

std::string PrimaryGeneratorAction::MakeCurrentSourceEventUid() const
{
  std::ostringstream oss;
  oss << MakeCurrentPhysicalEventUid() << "_"
      << fCurrentAlphaLiReplayIndex;
  return oss.str();
}

G4double PrimaryGeneratorAction::GetCurrentTrajectoryWeight() const
{
  return (fAlphaLiReplayPerCapture > 0)
             ? (1.0 / static_cast<G4double>(fAlphaLiReplayPerCapture))
             : 1.0;
}

std::string PrimaryGeneratorAction::MakeCurrentPhysicalEventUid() const
{
  std::ostringstream oss;
  oss << fCurrentRecord.eventID << "_"
      << fCurrentRecord.record_index;
  return oss.str();
}

// --------------------------------------------------------------------

std::string PrimaryGeneratorAction::DetermineSurfaceMode(const CaptureRecord &rec) const
{
  const auto *det = dynamic_cast<const DetectorConstruction *>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  const G4double screenT = rec.thickness_um;
  const G4double localT = det->GetEffectiveLocalThickness() / um;
  const G4double surfaceCut = 0.5 * localT;

  if (rec.depth_um <= surfaceCut)
  {
    return "front_surface";
  }

  if ((screenT - rec.depth_um) <= surfaceCut)
  {
    return "back_surface";
  }

  return "bulk";
}

// --------------------------------------------------------------------

G4double PrimaryGeneratorAction::DetermineTargetLocalZ(
    const CaptureRecord &rec, const std::string &mode) const
{
  const auto *det = dynamic_cast<const DetectorConstruction *>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  const G4double localT = det->GetEffectiveLocalThickness();

  if (mode == "front_surface")
  {
    return +0.5 * localT - rec.depth_um * um;
  }

  if (mode == "back_surface")
  {
    return -0.5 * localT + (rec.thickness_um - rec.depth_um) * um;
  }

  // Bulk events keep depth_um only as macroscopic metadata.
  // The local capture point is sampled later from the RVE bulk-safe BN volume.
  return 0.0;
}

// --------------------------------------------------------------------

G4bool PrimaryGeneratorAction::SelectBNSphereForTargetZ(
    G4double targetZ,
    G4ThreeVector &chosenCenter,
    G4double &usedZ,
    G4bool &usedFallback) const
{
  const auto *det = dynamic_cast<const DetectorConstruction *>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  const G4double R = det->GetBNRadius();

  auto trySelect = [&](const std::vector<G4ThreeVector> &centers,
                       G4bool requireSafeXY,
                       G4ThreeVector &centerOut,
                       G4double &zOut,
                       G4bool &fallbackOut) -> G4bool
  {
    std::vector<G4int> candidateIdx;
    std::vector<G4double> weights;

    for (G4int i = 0; i < static_cast<G4int>(centers.size()); ++i)
    {
      const G4double dz = targetZ - centers[i].z();
      if (std::abs(dz) < R)
      {
        const G4double baseW = pi * (R * R - dz * dz);
        if (baseW > 0.0)
        {
          const G4double diskR = std::sqrt(std::max(0.0, R * R - dz * dz));
          if (requireSafeXY && !IsSurfaceSliceSafelyInsideXY(centers[i], diskR, det))
          {
            continue;
          }

          if (!requireSafeXY && !DoesSurfaceSliceIntersectSafeXY(centers[i], diskR, det))
          {
            continue;
          }

          candidateIdx.push_back(i);
          weights.push_back(baseW);
        }
      }
    }

    if (!candidateIdx.empty())
    {
      const G4double sumW = std::accumulate(weights.begin(), weights.end(), 0.0);
      G4double pick = G4UniformRand() * sumW;

      for (G4int k = 0; k < static_cast<G4int>(candidateIdx.size()); ++k)
      {
        pick -= weights[k];
        if (pick <= 0.0 || k == static_cast<G4int>(candidateIdx.size()) - 1)
        {
          centerOut = centers[candidateIdx[k]];
          zOut = targetZ;
          fallbackOut = !requireSafeXY;
          return true;
        }
      }
    }

    return false;
  };

  // Prefer cross-section disks that remain well inside the RVE XY boundary.
  if (trySelect(det->GetSafeBNCenters(), true, chosenCenter, usedZ, usedFallback))
  {
    return true;
  }

  if (trySelect(det->GetPlacedBNCenters(), true, chosenCenter, usedZ, usedFallback))
  {
    return true;
  }

  // Last resort: preserve the surface-depth mapping, but only if the slice
  // intersects the XY safety window and a safe point can be rejection-sampled.
  if (trySelect(det->GetSafeBNCenters(), false, chosenCenter, usedZ, usedFallback))
  {
    return true;
  }

  if (trySelect(det->GetPlacedBNCenters(), false, chosenCenter, usedZ, usedFallback))
  {
    return true;
  }

  return false;
}

// --------------------------------------------------------------------

G4ThreeVector PrimaryGeneratorAction::SamplePointInSphereSlice(
    const G4ThreeVector &center,
    G4double zSlice,
    G4double sphereRadius) const
{
  const G4double dz = zSlice - center.z();
  const G4double diskR = std::sqrt(std::max(0.0, sphereRadius * sphereRadius - dz * dz));

  const G4double r = diskR * std::sqrt(G4UniformRand());
  const G4double phi = twopi * G4UniformRand();

  const G4double x = center.x() + r * std::cos(phi);
  const G4double y = center.y() + r * std::sin(phi);

  return G4ThreeVector(x, y, zSlice);
}

// --------------------------------------------------------------------

G4bool PrimaryGeneratorAction::SampleSafePointInSphereSlice(
    const G4ThreeVector &center,
    G4double zSlice,
    G4double sphereRadius,
    G4ThreeVector &point) const
{
  const auto *det = dynamic_cast<const DetectorConstruction *>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  for (G4int trial = 0; trial < 4096; ++trial)
  {
    const G4ThreeVector p = SamplePointInSphereSlice(center, zSlice, sphereRadius);
    if (IsInsideSafeXY(p, det))
    {
      point = p;
      return true;
    }
  }

  return false;
}

// --------------------------------------------------------------------

G4ThreeVector PrimaryGeneratorAction::SamplePointInSphereVolume(
    const G4ThreeVector &center,
    G4double sphereRadius) const
{
  const G4double cosTheta = 2.0 * G4UniformRand() - 1.0;
  const G4double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
  const G4double phi = twopi * G4UniformRand();
  const G4double r = sphereRadius * std::cbrt(G4UniformRand());

  return center + G4ThreeVector(
                      r * sinTheta * std::cos(phi),
                      r * sinTheta * std::sin(phi),
                      r * cosTheta);
}

// --------------------------------------------------------------------

G4bool PrimaryGeneratorAction::SampleBulkCapturePoint(
    G4ThreeVector &chosenCenter,
    G4ThreeVector &capturePoint,
    G4bool &usedFallback) const
{
  const auto *det = dynamic_cast<const DetectorConstruction *>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  const G4double R = det->GetBNRadius();

  auto collectFullySafeCenters = [&](const std::vector<G4ThreeVector> &centers)
  {
    std::vector<G4ThreeVector> out;
    for (const auto &c : centers)
    {
      if (IsBulkSphereFullyInSafeWindow(c, R, det))
      {
        out.push_back(c);
      }
    }
    return out;
  };

  auto chooseCenterUniformly = [](const std::vector<G4ThreeVector> &centers)
  {
    const G4int idx = std::min(
        static_cast<G4int>(centers.size()) - 1,
        static_cast<G4int>(G4UniformRand() * centers.size()));
    return centers[idx];
  };

  auto safeCenters = collectFullySafeCenters(det->GetSafeBNCenters());
  if (safeCenters.empty())
  {
    safeCenters = collectFullySafeCenters(det->GetPlacedBNCenters());
  }

  if (!safeCenters.empty())
  {
    chosenCenter = chooseCenterUniformly(safeCenters);
    capturePoint = SamplePointInSphereVolume(chosenCenter, R);
    usedFallback = false;
    return true;
  }

  // Fallback for sparse/pathological placements: find a BN sphere that still
  // contains at least one point in the safe bulk window, then rejection-sample.
  std::vector<G4ThreeVector> candidateCenters;
  for (const auto &c : det->GetSafeBNCenters())
  {
    if (DoesSphereIntersectBulkSafeWindow(c, R, det))
    {
      candidateCenters.push_back(c);
    }
  }

  if (candidateCenters.empty())
  {
    for (const auto &c : det->GetPlacedBNCenters())
    {
      if (DoesSphereIntersectBulkSafeWindow(c, R, det))
      {
        candidateCenters.push_back(c);
      }
    }
  }

  if (candidateCenters.empty())
  {
    return false;
  }

  chosenCenter = chooseCenterUniformly(candidateCenters);
  for (G4int trial = 0; trial < 4096; ++trial)
  {
    const G4ThreeVector p = SamplePointInSphereVolume(chosenCenter, R);
    if (IsBulkPointSafe(p, det))
    {
      capturePoint = p;
      usedFallback = true;
      return true;
    }
  }

  return false;
}

// --------------------------------------------------------------------

void PrimaryGeneratorAction::GenerateReactionProducts(
    G4Event *event,
    const G4ThreeVector &position,
    G4bool useGroundStateBranch) const
{
  auto *particleTable = G4ParticleTable::GetParticleTable();
  auto *ionTable = G4IonTable::GetIonTable();

  auto *alphaDef = particleTable->FindParticle("alpha");
  auto *li7Def = ionTable->GetIon(3, 7, 0.0);

  if (!alphaDef || !li7Def)
  {
    G4Exception("PrimaryGeneratorAction::GenerateReactionProducts",
                "BNZS007", FatalException,
                "Failed to obtain alpha or Li7 particle definition.");
    return;
  }

  G4double eAlpha = 0.0;
  G4double eLi7 = 0.0;

  // 10B(n,alpha)7Li
  // excited branch ~93.7%: alpha 1.47 MeV, Li7 0.84 MeV (+478 keV gamma, not emitted here)
  // ground branch  ~6.3% : alpha 1.78 MeV, Li7 1.01 MeV
  if (useGroundStateBranch)
  {
    eAlpha = 1.776 * MeV;
    eLi7 = 1.013 * MeV;
  }
  else
  {
    eAlpha = 1.470 * MeV;
    eLi7 = 0.840 * MeV;
  }

  const G4ThreeVector dir = RandomUnitVector();

  // alpha
  fParticleGun->SetParticleDefinition(alphaDef);
  fParticleGun->SetParticleEnergy(eAlpha);
  fParticleGun->SetParticlePosition(position);
  fParticleGun->SetParticleMomentumDirection(dir);
  fParticleGun->GeneratePrimaryVertex(event);

  // Li7 in opposite direction
  fParticleGun->SetParticleDefinition(li7Def);
  fParticleGun->SetParticleEnergy(eLi7);
  fParticleGun->SetParticlePosition(position);
  fParticleGun->SetParticleMomentumDirection(-dir);
  fParticleGun->GeneratePrimaryVertex(event);
}

// --------------------------------------------------------------------

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *event)
{
  const G4int geantEventID = event->GetEventID();

  if (fRemainingReplaysForCurrentCapture <= 0)
  {
    if (!PrepareCurrentCaptureReplayState())
    {
      G4cout
          << "[PrimaryGeneratorAction] No more input capture records. "
          << "Aborting run at Geant4 event " << geantEventID << G4endl;

      G4RunManager::GetRunManager()->AbortRun(true);
      return;
    }
  }

  fCurrentAlphaLiReplayIndex =
      std::max(0, fAlphaLiReplayPerCapture - fRemainingReplaysForCurrentCapture);

  // choose branch
  const G4bool useGroundStateBranch = (G4UniformRand() < 0.063);

  GenerateReactionProducts(
      event,
      fCurrentLocalCapturePosition,
      useGroundStateBranch);

  --fRemainingReplaysForCurrentCapture;

  if (geantEventID < 5)
  {
    G4cout
        << "\n[PrimaryGeneratorAction] Event " << geantEventID
        << "\n  input eventID      = " << fCurrentRecord.eventID
        << "\n  record index       = " << fCurrentRecord.record_index
        << "\n  replay index       = " << fCurrentAlphaLiReplayIndex
        << "\n  replay count       = " << fAlphaLiReplayPerCapture
        << "\n  mode               = " << fCurrentSurfaceMode
        << "\n  target z           = " << fCurrentTargetLocalZ / um << " um"
        << "\n  used z             = " << fCurrentUsedLocalZ / um << " um"
        << "\n  selected BN center = (" << fCurrentSelectedBNCenter.x() / um
        << ", " << fCurrentSelectedBNCenter.y() / um
        << ", " << fCurrentSelectedBNCenter.z() / um << ") um"
        << "\n  capture point      = (" << fCurrentLocalCapturePosition.x() / um
        << ", " << fCurrentLocalCapturePosition.y() / um
        << ", " << fCurrentLocalCapturePosition.z() / um << ") um"
        << G4endl;
  }
}
