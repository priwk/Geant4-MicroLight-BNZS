#ifndef StageDOpticalStats_h
#define StageDOpticalStats_h 1

#include "G4ThreeVector.hh"
#include "globals.hh"

#include <array>
#include <string>

struct StageDPhotonLaunchRecord
{
  G4int geantEventID = -1;
  G4double wavelength_nm = 0.0;
  std::string source_mode;
  std::string source_phase;
  G4ThreeVector source_position;
  G4ThreeVector momentum_direction;
  G4ThreeVector polarization;
  G4double photon_weight = 1.0;
  G4bool is_continuation = false;
};

struct StageDReentryDiagnosticRecord
{
  G4int event_id = -1;
  G4int reentry_index = 0;
  std::string strategy;
  std::string fallback_level;
  std::string exit_phase;
  std::string entry_phase;

  G4ThreeVector old_dir;
  G4ThreeVector exit_point;
  G4ThreeVector entry_point;

  G4double particle_q_exit = -1.0;
  G4double particle_q_entry = -1.0;
  G4double particle_mu_exit = 0.0;
  G4double particle_mu_entry = 0.0;

  G4double matrix_clearance_exit_um = -1.0;
  G4double matrix_clearance_entry_um = -1.0;
  std::string matrix_nearest_phase_exit;
  std::string matrix_nearest_phase_entry;
  G4int matrix_clearance_bin_exit = -1;
  G4int matrix_clearance_bin_entry = -1;

  G4int trials = 0;
};

struct StageDReentryPortalSummary
{
  static constexpr std::size_t kFaceCount = 6;
  static constexpr std::size_t kBinCount = 4;

  G4int total_portal_count = 0;
  std::array<G4int, kFaceCount> portal_count_by_face{};
  std::array<G4int, kBinCount> portal_count_by_bin{};
};

struct StageDPhotonEventRecord
{
  static constexpr std::size_t kPhaseFunctionBins = 64;

  G4int photonID = -1;
  std::string ratio;
  std::string placement_file;
  std::string source_mode;
  std::string boundary_mode;
  std::string reentry_mode;
  std::string particle_reentry_mode;
  std::string matrix_reentry_mode;
  G4double wavelength_nm = 0.0;

  std::string source_phase;
  G4double source_x_um = 0.0;
  G4double source_y_um = 0.0;
  G4double source_z_um = 0.0;

  std::string final_status;
  G4bool absorbed = false;

  G4double total_path_length_um = 0.0;
  G4double path_length_bn_um = 0.0;
  G4double path_length_zns_um = 0.0;
  G4double path_length_matrix_um = 0.0;
  G4double path_length_world_um = 0.0;
  G4int num_steps = 0;

  G4int num_absorbed_total = 0;
  G4int num_absorbed_BN = 0;
  G4int num_absorbed_ZnS = 0;
  G4int num_absorbed_Matrix = 0;
  G4int num_absorbed_World = 0;

  G4int num_encounter_total = 0;
  G4int num_encounter_BN = 0;
  G4int num_encounter_ZnS = 0;
  G4int num_encounter_effective_total = 0;
  G4int num_encounter_effective_BN = 0;
  G4int num_encounter_effective_ZnS = 0;
  G4double sum_cos_theta_encounter = 0.0;
  G4double sum_cos_theta_encounter_BN = 0.0;
  G4double sum_cos_theta_encounter_ZnS = 0.0;
  G4double sum_cos_theta_encounter_effective = 0.0;
  G4double sum_cos_theta_encounter_effective_BN = 0.0;
  G4double sum_cos_theta_encounter_effective_ZnS = 0.0;
  G4double sum_one_minus_cos_theta_encounter = 0.0;
  G4double sum_one_minus_cos_theta_encounter_BN = 0.0;
  G4double sum_one_minus_cos_theta_encounter_ZnS = 0.0;
  G4double sum_one_minus_cos_theta_encounter_effective = 0.0;
  G4double sum_one_minus_cos_theta_encounter_effective_BN = 0.0;
  G4double sum_one_minus_cos_theta_encounter_effective_ZnS = 0.0;
  G4double sum_cos2_theta_encounter = 0.0;
  G4double sum_cos2_theta_encounter_BN = 0.0;
  G4double sum_cos2_theta_encounter_ZnS = 0.0;
  G4double sum_cos2_theta_encounter_effective = 0.0;
  G4double sum_cos2_theta_encounter_effective_BN = 0.0;
  G4double sum_cos2_theta_encounter_effective_ZnS = 0.0;

  G4int num_particle_scatter = 0;
  G4int num_particle_scatter_BN = 0;
  G4int num_particle_scatter_ZnS = 0;
  G4int num_real_scatter = 0;
  G4int num_bulk_scatter = 0;
  G4int num_boundary_scatter = 0;
  G4int num_boundary_scatter_BN = 0;
  G4int num_boundary_scatter_ZnS = 0;
  G4int num_material_boundary = 0;
  G4int num_reentry = 0;
  G4int num_reentry_BN = 0;
  G4int num_reentry_ZnS = 0;
  G4int num_reentry_matrix = 0;
  G4int num_reentry_particle_q_mu = 0;
  G4int num_reentry_matrix_clearance_portal = 0;
  G4int num_reentry_fallback_same_bin = 0;
  G4int num_reentry_fallback_adjacent_bin = 0;
  G4int num_reentry_fallback_any_bin = 0;
  G4int num_reentry_fallback_any_phase_same_bin = 0;
  G4int num_reentry_fallback_any_portal = 0;
  G4int num_reentry_random_matrix_debug = 0;
  G4int num_reentry_failed = 0;

  G4double sum_cos_theta = 0.0;
  G4double sum_cos_theta_particle = 0.0;
  G4double sum_cos_theta_particle_BN = 0.0;
  G4double sum_cos_theta_particle_ZnS = 0.0;
  G4double sum_cos_theta_bulk = 0.0;
  G4double sum_cos_theta_boundary = 0.0;
  G4double sum_cos_theta_boundary_BN = 0.0;
  G4double sum_cos_theta_boundary_ZnS = 0.0;
  G4double g1_encounter_for_this_photon = 0.0;
  G4double g2_encounter_for_this_photon = 0.0;
  G4double mu_s_prime_direct_encounter_per_um_for_this_photon = 0.0;
  G4double g1_encounter_raw_for_this_photon = 0.0;
  G4double g2_encounter_raw_for_this_photon = 0.0;
  G4double mu_s_prime_direct_encounter_raw_per_um_for_this_photon = 0.0;
  G4double mean_cos_theta_particle_for_this_photon = 0.0;
  G4double mean_cos_theta_for_this_photon = 0.0;
  G4double mean_cos_theta_bulk_for_this_photon = 0.0;
  G4double mean_cos_theta_boundary_for_this_photon = 0.0;
  G4double weight = 1.0;

  std::array<G4int, kPhaseFunctionBins> phase_function_histogram{};

  G4bool in_particle_segment = false;
  std::string particle_segment_phase;
  G4ThreeVector particle_segment_entry_direction;
};

#endif
