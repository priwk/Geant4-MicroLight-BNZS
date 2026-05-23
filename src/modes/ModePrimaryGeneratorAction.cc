#include "ModePrimaryGeneratorAction.hh"

#include "AnalysisConfig.hh"
#include "PrimaryGeneratorAction.hh"
#include "StageAPrimaryGeneratorAction.hh"
#include "StageCOpticalPrimaryGeneratorAction.hh"
#include "StageDOpticalPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4Exception.hh"
#include "G4ios.hh"

ModePrimaryGeneratorAction::ModePrimaryGeneratorAction(AnalysisConfig *config)
    : G4VUserPrimaryGeneratorAction(),
      fConfig(config),
      fStageBPrimary(nullptr),
      fStageAPrimary(nullptr),
      fStageCPrimary(nullptr),
      fStageDPrimary(nullptr)
{
    if (fConfig == nullptr)
    {
        G4Exception("ModePrimaryGeneratorAction::ModePrimaryGeneratorAction",
                    "BNZS_MODE_PRI_001", FatalException,
                    "AnalysisConfig pointer is null.");
        return;
    }

    G4cout << "[ModePrimaryGeneratorAction] Dispatcher initialized."
           << " current runMode = "
           << AnalysisConfig::RunModeName(fConfig->runMode)
           << G4endl;
}

ModePrimaryGeneratorAction::~ModePrimaryGeneratorAction()
{
    delete fStageAPrimary;
    delete fStageBPrimary;
    delete fStageCPrimary;
    delete fStageDPrimary;
}

void ModePrimaryGeneratorAction::GeneratePrimaries(G4Event *event)
{
    if (fConfig == nullptr)
    {
        G4Exception("ModePrimaryGeneratorAction::GeneratePrimaries",
                    "BNZS_MODE_PRI_002", FatalException,
                    "AnalysisConfig pointer is null.");
        return;
    }

    switch (fConfig->runMode)
    {
    case RunMode::StageA_NeutronPatch:
        if (fStageAPrimary == nullptr)
        {
            fStageAPrimary = new StageAPrimaryGeneratorAction(fConfig);
        }
        fStageAPrimary->GeneratePrimaries(event);
        return;

    case RunMode::StageB_ReplayAlphaLi:
        if (fStageBPrimary == nullptr)
        {
            fStageBPrimary = new PrimaryGeneratorAction(fConfig);
        }
        fStageBPrimary->GeneratePrimaries(event);
        return;

    case RunMode::StageC_OpticalStub:
        G4Exception("ModePrimaryGeneratorAction::GeneratePrimaries",
                    "BNZS_MODE_PRI_005", FatalException,
                    "RunMode StageC_OpticalStub is selected, but Stage C primary generator is not implemented yet.");
        return;

    case RunMode::StageC_OpticalRVE:
        if (fStageCPrimary == nullptr)
        {
            fStageCPrimary = new StageCOpticalPrimaryGeneratorAction(fConfig);
        }
        fStageCPrimary->GeneratePrimaries(event);
        return;

    case RunMode::StageD_OpticalHomogenization:
        if (fStageDPrimary == nullptr)
        {
            fStageDPrimary = new StageDOpticalPrimaryGeneratorAction(fConfig);
        }
        fStageDPrimary->GeneratePrimaries(event);
        return;

    default:
        G4Exception("ModePrimaryGeneratorAction::GeneratePrimaries",
                    "BNZS_MODE_PRI_006", FatalException,
                    "Unknown run mode.");
        return;
    }
}

PrimaryGeneratorAction *ModePrimaryGeneratorAction::GetStageBPrimaryAction() const
{
    if (fStageBPrimary == nullptr &&
        fConfig != nullptr &&
        fConfig->runMode == RunMode::StageB_ReplayAlphaLi)
    {
        const_cast<ModePrimaryGeneratorAction *>(this)->fStageBPrimary =
            new PrimaryGeneratorAction(fConfig);
    }
    return fStageBPrimary;
}

StageCOpticalPrimaryGeneratorAction *ModePrimaryGeneratorAction::GetStageCPrimaryAction()
{
    if (fStageCPrimary == nullptr &&
        fConfig != nullptr &&
        fConfig->runMode == RunMode::StageC_OpticalRVE)
    {
        fStageCPrimary = new StageCOpticalPrimaryGeneratorAction(fConfig);
    }
    return fStageCPrimary;
}

StageDOpticalPrimaryGeneratorAction *ModePrimaryGeneratorAction::GetStageDPrimaryAction()
{
    if (fStageDPrimary == nullptr &&
        fConfig != nullptr &&
        fConfig->runMode == RunMode::StageD_OpticalHomogenization)
    {
        fStageDPrimary = new StageDOpticalPrimaryGeneratorAction(fConfig);
    }
    return fStageDPrimary;
}
