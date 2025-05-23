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
// G4CellScoreComposer
//
// Class description:
// This class will be created for every cell standard
// scoring should be applied. It does the actual scoring.
// GetStandardCellScoreValues() delivers the struct
// G4CellScoreValues  does calculations based on the
// sums of scores and delivers the results in
// G4CellScoreValues.
//
// Author: Michael Dressel, CERN
// --------------------------------------------------------------------
#ifndef G4CellScoreComposer_hh
#define G4CellScoreComposer_hh 1

#include "G4CellScoreValues.hh"

class G4Step;

class G4CellScoreComposer
{
 public:
  // get values for estimators based on track length
  void EstimatorCalculation(const G4Step& step);

  // Called when a track enters the cell
  void TrackEnters();

  // Called if the cell population is increased
  void NewTrackPopedUp();

  // Called for every collision in the cell, taking weight of the colliding particle
  void SetCollisionWeight(G4double weight);

  // Set the importance of the cell
  void SetImportnace(G4double importance);

  // Return scores stored in G4CellScoreValues
  const G4CellScoreValues& GetStandardCellScoreValues() const;

 private:
  mutable G4CellScoreValues fSCScoreValues;
};

std::ostream& operator<<(std::ostream& out, const G4CellScoreComposer& ps);

#endif
