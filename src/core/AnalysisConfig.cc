#include "AnalysisConfig.hh"

#include <cstdlib>
#include <cctype>
#include <string>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>

namespace
{
  std::filesystem::path DetectProjectRoot()
  {
    namespace fs = std::filesystem;

    std::error_code ec;
    const fs::path cwd = fs::current_path(ec);
    if (ec)
      return fs::path(".");

    std::vector<fs::path> candidates;
    candidates.push_back(cwd);
    if (cwd.filename() == "build" && cwd.has_parent_path())
      candidates.push_back(cwd.parent_path());
    if (cwd.has_parent_path())
      candidates.push_back(cwd.parent_path());

    for (const auto &candidate : candidates)
    {
      if (candidate.empty())
        continue;
      if (fs::exists(candidate / "CMakeLists.txt", ec) && !ec)
        return candidate;
      ec.clear();
      if (fs::exists(candidate / "README.md", ec) && !ec)
        return candidate;
      ec.clear();
    }

    return (cwd.filename() == "build" && cwd.has_parent_path()) ? cwd.parent_path() : cwd;
  }

  std::string ToLowerCopy(std::string s)
  {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char ch)
                   { return static_cast<char>(std::tolower(ch)); });
    return s;
  }

  bool IsTruthyEnv(const char *value)
  {
    if (value == nullptr)
      return false;

    const std::string v = ToLowerCopy(value);
    return v == "1" || v == "true" || v == "yes" || v == "on";
  }

  bool TryParseRatioFolderName(const std::string &name, double &bnWt, double &znsWt)
  {
    const std::size_t dashPos = name.find('-');
    if (dashPos == std::string::npos)
      return false;

    try
    {
      bnWt = std::stod(name.substr(0, dashPos));
      znsWt = std::stod(name.substr(dashPos + 1));
    }
    catch (...)
    {
      return false;
    }

    return bnWt > 0.0 && znsWt > 0.0;
  }

  std::vector<std::string> SplitCsvLine(const std::string &line)
  {
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string field;
    while (std::getline(ss, field, ','))
    {
      fields.push_back(field);
    }
    return fields;
  }

  bool ReadFirstStageCSourceMetadata(const std::string &path,
                                     double &bnWt,
                                     double &znsWt,
                                     std::string &placementFile)
  {
    std::ifstream fin(path.c_str());
    if (!fin)
      return false;

    std::string headerLine;
    if (!std::getline(fin, headerLine))
      return false;

    const auto headers = SplitCsvLine(headerLine);
    std::unordered_map<std::string, std::size_t> index;
    for (std::size_t i = 0; i < headers.size(); ++i)
    {
      index[headers[i]] = i;
    }

    const auto bnIt = index.find("bn_wt");
    const auto znsIt = index.find("zns_wt");
    const auto placementIt = index.find("placement_file");
    if (bnIt == index.end() || znsIt == index.end() || placementIt == index.end())
      return false;

    std::string line;
    while (std::getline(fin, line))
    {
      if (line.empty())
        continue;

      const auto fields = SplitCsvLine(line);
      if (fields.size() <= std::max({bnIt->second, znsIt->second, placementIt->second}))
        continue;

      try
      {
        bnWt = std::stod(fields[bnIt->second]);
        znsWt = std::stod(fields[znsIt->second]);
      }
      catch (...)
      {
        return false;
      }

      placementFile = fields[placementIt->second];
      return bnWt > 0.0 && znsWt > 0.0 && !placementFile.empty();
    }

    return false;
  }
}

std::filesystem::path AnalysisConfig::ProjectRootPath()
{
  return DetectProjectRoot();
}

std::string AnalysisConfig::PathForRecord(const std::filesystem::path &path)
{
  namespace fs = std::filesystem;

  if (path.empty())
    return "";

  const fs::path projectRoot = ProjectRootPath();
  std::error_code ec;

  fs::path normalized = path;
  if (path.is_absolute())
    normalized = fs::weakly_canonical(path, ec);
  else
    normalized = (projectRoot / path).lexically_normal();

  if (ec)
  {
    ec.clear();
    normalized = path.lexically_normal();
  }

  const fs::path normalizedRoot = fs::weakly_canonical(projectRoot, ec);
  if (ec)
  {
    ec.clear();
    return normalized.generic_string();
  }

  const fs::path relative = normalized.lexically_relative(normalizedRoot);
  if (!relative.empty() && relative.native().find("..") != 0)
    return relative.generic_string();

  return normalized.generic_string();
}

std::string AnalysisConfig::PathForRecord(const std::string &path)
{
  if (path.empty())
    return "";
  return PathForRecord(std::filesystem::path(path));
}

AnalysisConfig::AnalysisConfig()
    : runMode(RunMode::StageB_ReplayAlphaLi),
      patchXY_um(50.0),
      microThickness_um(30.0),
      bnWt(1.0),
      znsWt(2.0),
      useRandomPlacement(true),
      placementFilePath(""),
      captureCsvPath(""),
      captureInputDir(""),
      opticalSourcePath(""),
      sourceSampling("uniformAlongStep"),
      opticalSamplesPerStep(1),
      writeStageCPhotonCsv(false),
      opticalParamsProvided(false),
      opticalMatrixRIndex(0.0),
      opticalMatrixAbsLengthUm(5.0e5),
      opticalBnRIndex(2.1),
      opticalBnAbsLengthUm(3.0e4),
      opticalZnsRIndex(2.36),
      opticalZnsAbsLengthUm(1.5e3),
      stageD_wavelength_nm(450.0),
      stageD_source_mode("uniform_ZnS"),
      stageD_boundary_mode("same_phase_reentry"),
      stageD_reentry_mode("state_matched"),
      stageD_particle_reentry_mode("sphere_q_mu"),
      stageD_matrix_reentry_mode("clearance_binned_portal"),
      stageD_scatter_metric("particle_encounter_no_threshold"),
      stageD_target_primary_scatter(160),
      stageD_theta_threshold_deg(0.10),
      stageD_max_reentry(10000),
      stageD_max_steps(100000),
      stageD_max_path_length_um(5000.0),
      stageD_output_dir(""),
      stageD_portal_nu(128),
      stageD_portal_nv(128),
      stageD_portal_margin_um(0.0),
      stageD_clearance_bin0_um(0.02),
      stageD_clearance_bin1_um(0.10),
      stageD_clearance_bin2_um(0.50),
      stageD_max_particle_reentry_trials(64),
      stageD_max_portal_fallback_level(4),
      allowThicknessEqualLocalPatch(true)
{
  const char *runModeEnv = std::getenv("BNZS_RUN_MODE");
  if (runModeEnv != nullptr)
  {
    const std::string mode = ToLowerCopy(runModeEnv);
    if (mode == "stagea" || mode == "a" || mode == "stagea_neutronpatch")
    {
      runMode = RunMode::StageA_NeutronPatch;
    }
    else if (mode == "stageb" || mode == "b" || mode == "stageb_replayalphali")
    {
      runMode = RunMode::StageB_ReplayAlphaLi;
    }
    else if (mode == "stagec" || mode == "c" || mode == "stagec_opticalstub")
    {
      runMode = RunMode::StageC_OpticalStub;
    }
    else if (mode == "stagec_opticalrve" || mode == "opticalrve")
    {
      runMode = RunMode::StageC_OpticalRVE;
    }
    else if (mode == "staged_opticalhomogenization" || mode == "staged" || mode == "d")
    {
      runMode = RunMode::StageD_OpticalHomogenization;
    }
  }

  const char *placementEnv = std::getenv("BNZS_PLACEMENT_FILE");
  if (placementEnv != nullptr && std::string(placementEnv).size() > 0)
  {
    placementFilePath = placementEnv;
    useRandomPlacement = false;
  }

  const char *randomPlacementEnv = std::getenv("BNZS_USE_RANDOM_PLACEMENT");
  if (randomPlacementEnv != nullptr)
  {
    useRandomPlacement = IsTruthyEnv(randomPlacementEnv);
  }

  const char *captureCsvEnv = std::getenv("BNZS_INPUT_CSV");
  if (captureCsvEnv != nullptr && std::string(captureCsvEnv).size() > 0)
  {
    captureCsvPath = captureCsvEnv;

    double dirBnWt = 0.0;
    double dirZnsWt = 0.0;
    const std::string dirName = std::filesystem::path(captureCsvPath).parent_path().filename().string();
    if (TryParseRatioFolderName(dirName, dirBnWt, dirZnsWt))
    {
      bnWt = dirBnWt;
      znsWt = dirZnsWt;
    }
  }

  const char *captureDirEnv = std::getenv("BNZS_INPUT_DIR");
  if (captureDirEnv != nullptr && std::string(captureDirEnv).size() > 0)
  {
    captureInputDir = captureDirEnv;

    double dirBnWt = 0.0;
    double dirZnsWt = 0.0;
    const std::string dirName = std::filesystem::path(captureInputDir).filename().string();
    if (TryParseRatioFolderName(dirName, dirBnWt, dirZnsWt))
    {
      bnWt = dirBnWt;
      znsWt = dirZnsWt;
    }
  }

  const char *opticalSourceEnv = std::getenv("BNZS_OPTICAL_SOURCE_CSV");
  if (opticalSourceEnv == nullptr || std::string(opticalSourceEnv).empty())
  {
    opticalSourceEnv = std::getenv("BNZS_STAGEC_SOURCE_CSV");
  }
  if (opticalSourceEnv != nullptr && std::string(opticalSourceEnv).size() > 0)
  {
    opticalSourcePath = opticalSourceEnv;

    double sourceBnWt = 0.0;
    double sourceZnsWt = 0.0;
    std::string sourcePlacementFile;
    if (ReadFirstStageCSourceMetadata(opticalSourcePath,
                                      sourceBnWt,
                                      sourceZnsWt,
                                      sourcePlacementFile))
    {
      bnWt = sourceBnWt;
      znsWt = sourceZnsWt;
      placementFilePath = sourcePlacementFile;
      useRandomPlacement = false;
    }
  }

  const char *sourceSamplingEnv = std::getenv("BNZS_SOURCE_SAMPLING");
  if (sourceSamplingEnv != nullptr && std::string(sourceSamplingEnv).size() > 0)
  {
    sourceSampling = sourceSamplingEnv;
  }

  const char *samplesEnv = std::getenv("BNZS_SAMPLE_PHOTONS_PER_STEP");
  if (samplesEnv != nullptr && std::string(samplesEnv).size() > 0)
  {
    try
    {
      const int samples = std::stoi(samplesEnv);
      if (samples > 0)
      {
        opticalSamplesPerStep = samples;
      }
    }
    catch (...)
    {
    }
  }

  const char *writePhotonCsvEnv = std::getenv("BNZS_WRITE_STAGEC_PHOTON_CSV");
  if (writePhotonCsvEnv != nullptr)
  {
    writeStageCPhotonCsv = IsTruthyEnv(writePhotonCsvEnv);
  }

  auto readOpticalDoubleEnv = [&](const char *name, double &target)
  {
    const char *value = std::getenv(name);
    if (value == nullptr || std::string(value).empty())
      return;
    try
    {
      const double parsed = std::stod(value);
      if (parsed > 0.0)
      {
        target = parsed;
        opticalParamsProvided = true;
      }
    }
    catch (...)
    {
    }
  };

  readOpticalDoubleEnv("BNZS_OPTICAL_MATRIX_RINDEX", opticalMatrixRIndex);
  readOpticalDoubleEnv("BNZS_OPTICAL_MATRIX_ABSLENGTH_UM", opticalMatrixAbsLengthUm);
  readOpticalDoubleEnv("BNZS_OPTICAL_BN_RINDEX", opticalBnRIndex);
  readOpticalDoubleEnv("BNZS_OPTICAL_BN_ABSLENGTH_UM", opticalBnAbsLengthUm);
  readOpticalDoubleEnv("BNZS_OPTICAL_ZNS_RINDEX", opticalZnsRIndex);
  readOpticalDoubleEnv("BNZS_OPTICAL_ZNS_ABSLENGTH_UM", opticalZnsAbsLengthUm);

  readOpticalDoubleEnv("BNZS_STAGED_WAVELENGTH_NM", stageD_wavelength_nm);
  readOpticalDoubleEnv("BNZS_STAGED_THETA_THRESHOLD_DEG", stageD_theta_threshold_deg);
  readOpticalDoubleEnv("BNZS_STAGED_MAX_PATH_LENGTH_UM", stageD_max_path_length_um);

  auto readPositiveIntEnv = [](const char *name, int &target)
  {
    const char *value = std::getenv(name);
    if (value == nullptr || std::string(value).empty())
      return;
    try
    {
      const int parsed = std::stoi(value);
      if (parsed > 0)
        target = parsed;
    }
    catch (...)
    {
    }
  };

  readPositiveIntEnv("BNZS_STAGED_MAX_REENTRY", stageD_max_reentry);
  readPositiveIntEnv("BNZS_STAGED_MAX_STEPS", stageD_max_steps);
  readPositiveIntEnv("BNZS_STAGED_TARGET_PRIMARY_SCATTER", stageD_target_primary_scatter);
  readPositiveIntEnv("BNZS_STAGED_PORTAL_NU", stageD_portal_nu);
  readPositiveIntEnv("BNZS_STAGED_PORTAL_NV", stageD_portal_nv);
  readPositiveIntEnv("BNZS_STAGED_MAX_PARTICLE_REENTRY_TRIALS", stageD_max_particle_reentry_trials);
  readPositiveIntEnv("BNZS_STAGED_MAX_PORTAL_FALLBACK_LEVEL", stageD_max_portal_fallback_level);

  const char *stageDSourceModeEnv = std::getenv("BNZS_STAGED_SOURCE_MODE");
  if (stageDSourceModeEnv != nullptr && std::string(stageDSourceModeEnv).size() > 0)
    stageD_source_mode = stageDSourceModeEnv;

  const char *stageDBoundaryModeEnv = std::getenv("BNZS_STAGED_BOUNDARY_MODE");
  if (stageDBoundaryModeEnv != nullptr && std::string(stageDBoundaryModeEnv).size() > 0)
    stageD_boundary_mode = stageDBoundaryModeEnv;

  const char *stageDReentryModeEnv = std::getenv("BNZS_STAGED_REENTRY_MODE");
  if (stageDReentryModeEnv != nullptr && std::string(stageDReentryModeEnv).size() > 0)
    stageD_reentry_mode = stageDReentryModeEnv;

  const char *stageDParticleReentryModeEnv = std::getenv("BNZS_STAGED_PARTICLE_REENTRY_MODE");
  if (stageDParticleReentryModeEnv != nullptr && std::string(stageDParticleReentryModeEnv).size() > 0)
    stageD_particle_reentry_mode = stageDParticleReentryModeEnv;

  const char *stageDMatrixReentryModeEnv = std::getenv("BNZS_STAGED_MATRIX_REENTRY_MODE");
  if (stageDMatrixReentryModeEnv != nullptr && std::string(stageDMatrixReentryModeEnv).size() > 0)
    stageD_matrix_reentry_mode = stageDMatrixReentryModeEnv;

  const char *stageDScatterMetricEnv = std::getenv("BNZS_STAGED_SCATTER_METRIC");
  if (stageDScatterMetricEnv != nullptr && std::string(stageDScatterMetricEnv).size() > 0)
    stageD_scatter_metric = stageDScatterMetricEnv;

  const char *stageDOutputDirEnv = std::getenv("BNZS_STAGED_OUTPUT_DIR");
  if (stageDOutputDirEnv != nullptr && std::string(stageDOutputDirEnv).size() > 0)
    stageD_output_dir = stageDOutputDirEnv;

  readOpticalDoubleEnv("BNZS_STAGED_PORTAL_MARGIN_UM", stageD_portal_margin_um);
  readOpticalDoubleEnv("BNZS_STAGED_CLEARANCE_BIN0_UM", stageD_clearance_bin0_um);
  readOpticalDoubleEnv("BNZS_STAGED_CLEARANCE_BIN1_UM", stageD_clearance_bin1_um);
  readOpticalDoubleEnv("BNZS_STAGED_CLEARANCE_BIN2_UM", stageD_clearance_bin2_um);

  const char *bnWtEnv = std::getenv("BNZS_BN_WT");
  const char *znsWtEnv = std::getenv("BNZS_ZNS_WT");
  if (bnWtEnv != nullptr && znsWtEnv != nullptr)
  {
    try
    {
      const double envBnWt = std::stod(bnWtEnv);
      const double envZnsWt = std::stod(znsWtEnv);
      if (envBnWt > 0.0 && envZnsWt > 0.0)
      {
        bnWt = envBnWt;
        znsWt = envZnsWt;
      }
    }
    catch (...)
    {
    }
  }
}

AnalysisConfig::~AnalysisConfig() = default;

const char *AnalysisConfig::RunModeName(RunMode mode)
{
  switch (mode)
  {
  case RunMode::StageA_NeutronPatch:
    return "StageA_NeutronPatch";
  case RunMode::StageB_ReplayAlphaLi:
    return "StageB_ReplayAlphaLi";
  case RunMode::StageC_OpticalStub:
    return "StageC_OpticalStub";
  case RunMode::StageC_OpticalRVE:
    return "StageC_OpticalRVE";
  case RunMode::StageD_OpticalHomogenization:
    return "StageD_OpticalHomogenization";
  default:
    return "UnknownRunMode";
  }
}
