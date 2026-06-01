#include "PhysicsList.hh"

#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4SystemOfUnits.hh"

PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
  SetVerboseLevel(1);

  // 1. 电磁过程：高精度低能带电粒子（alpha / Li ion）输运与能量沉积
  RegisterPhysics(new G4EmStandardPhysics_option4());

  // Optical physics stays globally available because Stage C/D select their
  // modes via macro commands before /run/initialize, after the physics list
  // object is already constructed.
  RegisterPhysics(new G4OpticalPhysics());

  // 3. 强子过程：高精度热中子物理（包含中子弹性散射、非弹性散射、捕获和裂变）
  RegisterPhysics(new G4HadronElasticPhysicsHP());
  RegisterPhysics(new G4HadronPhysicsQGSP_BIC_HP());

  // 4. 离子与停止过程：确保次级反冲核（如 Li 离子）的正确输运和停止行为
  RegisterPhysics(new G4IonPhysics());
  RegisterPhysics(new G4StoppingPhysics());

  // 5. 衰变过程（可选，看你的模拟是否涉及放射性同位素衰变）
  RegisterPhysics(new G4DecayPhysics());
}

PhysicsList::~PhysicsList() {}

void PhysicsList::SetCuts()
{
  // 阈值设置：0.1 um 对于微米级发光材料非常合适，能精确追踪 alpha 粒子
  SetCutValue(1.0 * um, "gamma");
  SetCutValue(0.1 * um, "e-");
  SetCutValue(0.1 * um, "e+");
  SetCutValue(0.1 * um, "proton");

  G4VUserPhysicsList::SetCuts();
}
