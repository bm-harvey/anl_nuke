// -*- mode: c++ -*-
//
// $Id: FaustFilter.h,v 1.2 2022/11/07 17:45:50 andyhannaman Exp $
// $Author: andyhannaman $
// $Date: 2022/11/07 17:45:50 $
//
#ifndef E101014_FaustFilter
#define E101014_FaustFilter
#ifndef FAUST_FaustParameterManager
#include "FaustParameterManager.h"
#endif

#include "TObjArray.h"
#include "BrPlane3D.h"

class FaustFilterParticle : public TObject
{
public:
   FaustFilterParticle();
   ~FaustFilterParticle() {}
   void Reset();
   void AddHit(Int_t idet, TObject *partObj);

   // Setters
   void SetPunchesThrough(Bool_t above = kTRUE) { fPunchesThrough = above; }
   void SetAboveThreshUnderSatReject(Bool_t reject = kTRUE) { fAboveThreshUnderSatReject = reject; }
   void SetFilteredE(Double_t e) { fFilteredE = e; }
   void SetFilteredTheta(Double_t theta) { fFilteredTheta = theta; }
   void SetFilteredPhi(Double_t phi) { fFilteredPhi = phi; }
   void SetFilteredXPos(Double_t xpos) { fXPos = xpos; }
   void SetFilteredYPos(Double_t ypos) { fYPos = ypos; }
   void SetSiEnergy(Double_t SiE) { fSiE = SiE; }
   void SetCsIEnergy(Double_t CsIE) { fCsIE = CsIE; }
   void SetXPos(Double_t XPos) { fXPos_U = XPos; }
   void SetYPos(Double_t YPos) { fYPos_U = YPos; }
   void SetTheta(Double_t Theta) { fTheta = Theta; }
   void SetPhi(Double_t Phi) { fPhi = Phi; }

   // Getters
   TObject *GetObjectAt(Int_t i);
   Bool_t ValidHit();
   Int_t GetIDet() const { return fiDet; }
   Int_t GetNumHits() const { return fNumHits; }
   Double_t GetFilteredE() const { return fFilteredE; }
   Double_t GetFilteredTheta() const { return fFilteredTheta; }
   Double_t GetFilteredPhi() const { return fFilteredPhi; }
   Double_t GetSiEnergy() const { return fSiE; }
   Double_t GetCsIEnergy() const { return fCsIE; }
   Double_t GetFilteredXPos() const {return fXPos;}
   Double_t GetFilteredYPos() const {return fYPos;}
   Double_t GetXPos() const {return fXPos_U;}
   Double_t GetYPos() const {return fYPos_U;}
   Double_t GetTheta() const {return fTheta;}
   Double_t GetYPhi() const {return fPhi;}
   

private:
   
   Int_t fiDet;
   Double_t fTheta;
   Double_t fPhi;
   Double_t fXPos_U;
   Double_t fYPos_U;
   Double_t fFilteredE;
   Double_t fFilteredTheta;
   Double_t fFilteredPhi;
   Double_t fXPos;
   Double_t fYPos;
   Double_t fSiE;
   Double_t fCsIE;
   Bool_t fPunchesThrough;
   Bool_t fAboveThreshUnderSatReject; // if true, rejected for AboveThreshUnderSat
   Int_t fNumHits;
   TObjArray fObjectList;
   TObjArray fHitList;
   
};

class CycSrimDev;
class CycDetectorVolume;

class FaustFilter : public TObject
{
public:
   // singleton
   static FaustFilter *Instance();
   virtual ~FaustFilter();

   void Init();
   void Begin();
   void Begin(Int_t run);
   void Clear(Option_t *option = "");

   Int_t Filter(Int_t z, Int_t a, Double_t e, Double_t theta, Double_t phi, TObject *partObj = 0);

   Int_t GeoFilter(Double_t Q, Double_t F); // theta and phi in DEGREES
   Int_t GeoFilter(Double_t Q, Double_t F, Double_t &intgerX, Double_t &intgerY, Double_t &intgerZ);
   void InitGeometry();
   void ReadGeometry();
   void CalcPlanes();
   TObjArray *GetListOfHits();
   TObjArray *GetListOfValidHits();
   TObjArray *GetListOfValidFilteredParticles();


   void IncludePositionResolution(Bool_t enable = kTRUE) { fIncludePositionResolution = enable; }
   void IncludeThresholdSaturation(Bool_t enable = kTRUE) { fIncludeThresholdSaturation = enable; }
   void IncludeEnergySmear(Bool_t enable = kTRUE) { fIncludeEnergySmear = enable; }
   void IncludePunchthrough(Bool_t enable = kTRUE) { fIncludePunchthrough = enable; }
   void IncludeGeometryOnly(Bool_t enable = kFALSE) { fIncludeGeometryOnly = enable; }
   void IncludeDetectorPositionSmear(Bool_t enable = kFALSE) { fIncludeDetectorPositionSmear = enable; }
   void SetDetPosTranslationSmear(Double_t smear = 0.1) { fTranslationSmear = smear; } //set the parameter for how much +/- to smear in x and y for each det position (cm)
   void SetDetPosRotationSmear(Double_t smear = 5.) { fRotationSmear = smear; } //set the parameter for how much +/- to smear rotation for each det position (degrees)
   void SetSiSmear(Double_t smear = 1.) { fSiSmear = smear; } //percent resolution for Si smearing (FWHM)
   void SetCsISmear(Double_t smear = 1.5) { fCsISmear = smear; } //percent resolution for CsISmearing (FWHM)


private:
   FaustFilter();

   static FaustFilter *fgInstance;
   Bool_t fIsInitialized;

   Double_t fCorner[69][4][3];
   Double_t fPlane[69][4];
   Double_t fEnergyLoss;
   FaustParameterManager *fParameterManager;

   void CalcInter(Double_t Q, Double_t F, Int_t det, Double_t &X, Double_t &Y, Double_t &Z);
   Double_t AreaTri(Double_t x0, Double_t y0, Double_t z0, Int_t det, Int_t c1, Int_t c2);

   Int_t IsHit(Double_t theta, Double_t phi, Double_t &posX, Double_t &posY, Double_t &posZ, Double_t &interX, Double_t &interY, Double_t &interZ);
   Bool_t PunchesThrough(Int_t idet, Int_t z, Int_t a, Double_t energy);
   Bool_t AboveThreshUnderSat(Int_t detId, Int_t z, Int_t a, Double_t energy, Double_t posX, Double_t posY);
   void SplitSignals(Double_t energyloss, Double_t XPos, Double_t YPos, Double_t &F1, Double_t &F2, Double_t &B1, Double_t &B2);
   void GetSmearedPosition(Int_t idet, Double_t energyloss, Double_t XPos, Double_t YPos, Double_t &XPos_S, Double_t &YPos_S);
   void GetSiEnergyLoss(Int_t idet, Double_t energy, Int_t z, Int_t a);
   void GetGlobalVector(Int_t idet, Double_t PosX, Double_t PosY, Double_t PosZ, Double_t &interX, Double_t &interY, Double_t &interZ);
   Double_t GetCsIPunchThrough(Int_t idet, Double_t z, Double_t a);
   Double_t GetSmearedSiEnergy();
   Double_t GetSmearedCsIEnergy(Double_t energy, Double_t MylarLoss);
   Double_t GetMylarEnergyLoss(Int_t idet, Double_t energy, Double_t z, Double_t a);
   void GetFaustThetaPhi(Double_t interX, Double_t interY, Double_t interZ, Double_t &theta, Double_t &phi);

   CycSrimDev *fSrim[68];
   CycSrimDev *fSrimCsI[68];
   CycSrimDev *fSrimMylar[68];
   CycDetectorVolume *fDetVol[68];
   BrPlane3D fFrontPlane[68];
   Double_t fSiPunchthrough[68][5];    // precalculated for p,d,t,3He,4He
   double_t fSiCsIPunchthrough[68][5]; // precalculated for pdt 3He 4He
   Double_t fXSatSlope[68];
   Double_t fXSatOffset[68];
   Double_t fYSatSlope[68];
   Double_t fYSatOffset[68];
   Double_t fThreshold[68];
   Double_t fTriggerNoise[68];
   Double_t fUncorrNoise[68];

   Int_t fRun;

   FaustFilterParticle fFilteredParticleList[68];
  Double_t fXDetSmear[68];
  Double_t fYDetSmear[68];
  Double_t fRotDetSmear[68];
   TObjArray fHitList; // used to return list of valid hits
   Bool_t fIncludePositionResolution;
   Bool_t fIncludeThresholdSaturation;
   Bool_t fIncludeEnergySmear;
   Bool_t fIncludePunchthrough;
   Bool_t fIncludeGeometryOnly;
   Bool_t fIncludeDetectorPositionSmear;
   Double_t fTranslationSmear;
   Double_t fRotationSmear;
   Double_t fSiSmear;
   Double_t fCsISmear;


public:
   ClassDef(FaustFilter, 0)

}; // end FaustFilter declaration

#endif
