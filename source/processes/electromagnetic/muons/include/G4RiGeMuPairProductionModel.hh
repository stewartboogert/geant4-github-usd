//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4RiGeMuPairProductionModel
//
// Authors:       Girardo Depaola & Ricardo Pacheco
//
// Creation date: 29.10.2024
//
//
// -------------------------------------------------------------------
//

#ifndef G4RiGeMuPairProductionModel_h
#define G4RiGeMuPairProductionModel_h 1

#include "G4VEmModel.hh"
#include "G4NistManager.hh"
#include "G4ElementData.hh"
#include "G4Physics2DVector.hh"
#include "G4VEmAngularDistribution.hh"
#include <vector>

class G4Element;
class G4ParticleChangeForLoss;
class G4RiGeAngularGenerator;

class G4RiGeMuPairProductionModel : public G4VEmModel
{
public:

  explicit G4RiGeMuPairProductionModel(const G4ParticleDefinition* p = nullptr);

  ~G4RiGeMuPairProductionModel() override = default;

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  void InitialiseLocal(const G4ParticleDefinition*,
                       G4VEmModel* masterModel) override;
      
  G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
              G4double kineticEnergy,
              G4double Z, G4double A,
              G4double cutEnergy,
              G4double maxEnergy) override;
         
  G4double ComputeDEDXPerVolume(const G4Material*,
                                const G4ParticleDefinition*,
                                G4double kineticEnergy,
                                G4double cutEnergy) override;

  void SampleSecondaries(std::vector<G4DynamicParticle*>*, 
                         const G4MaterialCutsCouple*,
                         const G4DynamicParticle*,
                         G4double tmin,
                         G4double maxEnergy) override;

  G4double MinPrimaryEnergy(const G4Material*,
                            const G4ParticleDefinition*,
                            G4double) override;

  G4double 
  ComputeDMicroscopicCrossSection(G4double tkin, G4double Z,
          G4double pairEnergy);

  inline void SetLowestKineticEnergy(G4double e);

  inline void SetParticle(const G4ParticleDefinition*);

  // hide assignment operator and copy constructor
  G4RiGeMuPairProductionModel& operator=
  (const G4RiGeMuPairProductionModel& right) = delete;
  G4RiGeMuPairProductionModel(const G4RiGeMuPairProductionModel&) = delete;

protected:

  G4double ComputMuPairLoss(G4double Z, G4double tkin, G4double cut,
                            G4double tmax);

  G4double ComputeMicroscopicCrossSection(G4double tkin,
                                          G4double Z,
                                          G4double cut);

  G4double FindScaledEnergy(G4int Z, G4double rand, G4double logTkin,
          G4double yymin, G4double yymax); 

  inline G4double MaxSecondaryEnergyForElement(G4double kineticEnergy,
                 G4double Z);

  void MakeSamplingTables();

  void StoreTables() const;

  G4bool RetrieveTables();

  virtual void DataCorrupted(G4int Z, G4double logTkin) const;

  G4ParticleChangeForLoss* fParticleChange = nullptr;
  const G4ParticleDefinition* particle = nullptr;
  G4NistManager* nist = nullptr;

  G4double factorForCross;
  G4double sqrte;
  G4double particleMass = 0.0;
  G4double z13 = 0.0;
  G4double z23 = 0.0;
  G4double lnZ = 0.0;

  G4double minPairEnergy;
  G4double lowestKinEnergy;

  G4double emin;
  G4double emax;
  G4double ymin = -5.0;
  G4double dy = 0.005;

  // Random numbers for sampling
  G4double randNumbs[9];

  G4int currentZ = 0;
  G4int nYBinPerDecade = 4;
  std::size_t nbiny = 1000;
  std::size_t nbine = 0;

  G4bool fTableToFile = false;

  // static members
  static const G4int NZDATPAIR = 5;
  static const G4int NINTPAIR = 8;
  static const G4int ZDATPAIR[NZDATPAIR];
  static const G4double xgi[NINTPAIR];
  static const G4double wgi[NINTPAIR];

private:

  G4RiGeAngularGenerator* fAngularGenerator;
  G4ParticleDefinition* theElectron;
  G4ParticleDefinition* thePositron;
  G4String dataName{""};
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4RiGeMuPairProductionModel::SetLowestKineticEnergy(G4double e) 
{
  lowestKinEnergy = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void G4RiGeMuPairProductionModel::SetParticle(const G4ParticleDefinition* p)
{
  if(nullptr == particle) {
    particle = p;
    particleMass = particle->GetPDGMass();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double 
G4RiGeMuPairProductionModel::MaxSecondaryEnergyForElement(G4double kineticEnergy,
                         G4double ZZ)
{
  G4int Z = G4lrint(ZZ);
  if(Z != currentZ) {
    currentZ = Z;
    z13 = nist->GetZ13(Z);
    z23 = z13*z13;
    lnZ = nist->GetLOGZ(Z);
  }
  return kineticEnergy + particleMass*(1.0 - 0.75*sqrte*z13);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
