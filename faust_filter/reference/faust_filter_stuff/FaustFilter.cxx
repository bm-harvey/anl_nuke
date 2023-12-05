// $Id: FaustFilter.cxx,v 1.2 2022/11/07 17:45:41 andyhannaman Exp $
// $Author: andyhannaman $
// $Date: 2022/11/07 17:45:41 $
//
#include "FaustFilter.h"

#include <math.h>

#include <fstream>
#include <iostream>

#include "CycDetectorVolume.h"
#include "CycSrimDev.h"
#include "PathManager.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TString.h"
using namespace std;

FaustFilterParticle::FaustFilterParticle()
{
  fObjectList.SetOwner(kFALSE);
  fPunchesThrough = kFALSE;
  fAboveThreshUnderSatReject = kFALSE;
}

void FaustFilterParticle::Reset()
{
  fNumHits = 0;
  fObjectList.Clear(); // should not delete since owner is kFALSE
  fFilteredE = 0.0;
  fFilteredTheta = 0.0;
  fFilteredPhi = 0.0;
  fPunchesThrough = kFALSE;
  fAboveThreshUnderSatReject = kFALSE; // assume most will make efficiency cut
}

void FaustFilterParticle::AddHit(Int_t idet, TObject *partObj)
{
  fNumHits++;
  fiDet = idet;
  if (partObj)
    fObjectList.Add(partObj);
}

TObject *FaustFilterParticle::GetObjectAt(Int_t i)
{
  return fObjectList.At(i);
}

Bool_t FaustFilterParticle::ValidHit()
{
  // if(fNumHits > 1) printf("ValidHit: detid = %d, num hits = %d, aboveThresh = %d, eff = %d\n",fDetId,fNumHits,fPunchesThrough,fEfficiencyReject);
  if (fNumHits == 1 && fPunchesThrough && !fAboveThreshUnderSatReject)
    return kTRUE;
  else
    return kFALSE;
}

//==============================================================================

FaustFilter *FaustFilter::fgInstance = 0;

FaustFilter::FaustFilter()
{
  InitGeometry();
  fRun = 0; // Not using for the moment
  fHitList.SetOwner(kFALSE);
  fHitList.Clear(); // start on right track

  fIncludePositionResolution = kTRUE;  // default
  fIncludeThresholdSaturation = kTRUE; // default
  fIncludeEnergySmear = kTRUE;         // default
  fIncludePunchthrough = kTRUE;        // default
  fIncludeGeometryOnly = kFALSE;               // default is FALSE
  fIncludeDetectorPositionSmear = kFALSE;               // default is FALSE
  fRotationSmear = 5.; //degrees
  fTranslationSmear = 0.1; //cm
  fSiSmear = 1.;
  fCsISmear = 1.5;
}

FaustFilter::~FaustFilter()
{
}

FaustFilter *FaustFilter::Instance()
{
  // Create an instance of this object

  if (!fgInstance)
    fgInstance = new FaustFilter();
  return fgInstance;
}

//==============================================================================
void FaustFilter::Init()
{

  gRandom = new TRandom3(0); // initialize TRandom3 used in threshold noise

  // start up the parameter manager
  fParameterManager = FaustParameterManager::Instance();
  fParameterManager->Init();
  // Setup srim for calculating thresholds
  // Calculating one for each detector in case we put the different
  // thicknesses for the different detectors
  const Int_t zpar[5] = {1, 1, 1, 2, 2};
  const Int_t apar[5] = {1, 2, 3, 3, 4};

  //set the x and y translation and detector rotation values for entire array - retain the detector plane orientation facing target. - this approximates the error in the dadl PCB screw holes.// for now 1 mm x and y and 5 degrees

    for (Int_t i = 0; i < 68; i++)
  {
    fXDetSmear[i] = ((gRandom->Rndm()*fTranslationSmear*2.) - fTranslationSmear);
    // cout <<"xsmear " << i << ": " << fXDetSmear[i] << endl;
    fYDetSmear[i] = ((gRandom->Rndm()*fTranslationSmear*2.) - fTranslationSmear);
    fRotDetSmear[i] = TMath::DegToRad()*((gRandom->Rndm()*fRotationSmear*2.) - fRotationSmear); //+/- 5 degrees
  }

  // add section to read in txt file with DADL threshold vs. x,y position and saturation vs. x,y, position

  // Get detector thickness from parameter manager

  for (Int_t i = 0; i < 68; i++)
  {

    Double_t thick = fParameterManager->GetDetThickness(i); // returns detector thickness in um
 
    CycSrimDev *srim;
    srim = new CycSrimDev(CycSrimDev::SrimMaterialSi, thick, CycSrimDev::kUnitsMicron); // need to use CycSrimDevdev
    fSrim[i] = srim;

    // calculate CsI punch through energy for the idet, z, and a
    Double_t CsILength = fParameterManager->GetCsILength(i); // cm
    CsILength = CsILength * 1e4;
    CycSrimDev *srimCsI;
    srimCsI = new CycSrimDev(CycSrimDev::SrimMaterialCsI, CsILength, CycSrimDev::kUnitsMicron);
    fSrimCsI[i] = srimCsI;

    Double_t mylarthick = fParameterManager->GetMylarThickness(i);
    CycSrimDev *srimMylar;
    srimMylar = new CycSrimDev(CycSrimDev::SrimMaterialMylar, mylarthick, CycSrimDev::kUnitsMgCm2);
    fSrimMylar[i] = srimMylar;

    TString outLine = Form("det %d", i);
    for (Int_t ipar = 0; ipar < 5; ipar++)
    {
      Int_t z = zpar[ipar];
      Int_t a = apar[ipar];
      fSiPunchthrough[i][ipar] = srim->GetPunchthroughEnergy(z, a);
      fSiCsIPunchthrough[i][ipar] = GetCsIPunchThrough(i, z, a);
      outLine += Form("  %f  ", fSiPunchthrough[i][ipar]);
    }
    // printf("%s\n",outLine.Data());
  }

  // read in signal thresholds  obtained using DADL box model (/data/sjygroup/sjy42/060821/HistogramScripts/MatchThreshold_Sig.C) and saturations obtained from exp data

  PathManager *pathMan = PathManager::Instance();
  const Char_t *dataDir = pathMan->GetDataDir();

  ifstream inFileThresh;
  ifstream inFileSat;
  // inFileThresh.open("/home/sjygroup/ahannaman/060821Calib/CalibrationParameters/FaustDetectorParameters/detnum_AMDThreshold.txt", ios::in);
  // inFileSat.open("/home/sjygroup/ahannaman/060821Calib/CalibrationParameters/FaustDetectorParameters/detnum_AMDSaturation.txt", ios::in);
  inFileThresh.open(Form("%s/faustGeom/detnum_AMDThreshold.txt", dataDir), ios::in);
  inFileSat.open(Form("%s/faustGeom/detnum_AMDSaturation.txt", dataDir), ios::in);
  Double_t Thresh = 0;
  Double_t Trig = 0;
  Double_t Uncorr = 0;
  Double_t XSlope = 0;
  Double_t YSlope = 0;
  Double_t XInt = 0;
  Double_t YInt = 0;
  TString lineThresh, lineSat, line;
  for (Int_t i = 0; i < 70; i++)
  {
    lineThresh.ReadLine(inFileThresh);
    lineSat.ReadLine(inFileSat);

    if (inFileThresh.eof() && inFileSat.eof())
      break;

    if (!lineThresh.BeginsWith("#"))
    {
      Int_t id = 397388833; // needs to be ridiculous number since det # starts from 0 and it may not match

      sscanf(lineThresh.Data(), "%d %lf %lf %lf", &id, &Thresh, &Trig, &Uncorr);
      // cout <<"id: " << id << endl;
      // cout <<"Thresh: " << Thresh << endl;
      fThreshold[id] = Thresh;
      fTriggerNoise[id] = Trig;
      fUncorrNoise[id] = Uncorr;
      // cout <<"id: " << id << endl;
      // cout <<"Thresh: " << Thresh << endl;
      // cout <<"Noise: " << Trig << endl;
      // cout <<"Unc: " << Uncorr << endl;
    }
    if (!lineSat.BeginsWith("#"))
    {
      Int_t id = 397388833; // needs to be ridiculous number since det # starts from 0 and it may not match

      sscanf(lineSat.Data(), "%d %lf %lf %lf %lf", &id, &XSlope, &XInt, &YSlope, &YInt);

      fXSatSlope[id] = XSlope;
      fXSatOffset[id] = XInt;
      fYSatSlope[id] = YSlope;
      fYSatOffset[id] = YInt;
      // cout <<"id: " << id << endl;
      // cout <<"XSl: " << XSlope << endl;
      // cout <<"YSl: " << YSlope << endl;
      // cout <<"YInt: " << YInt << endl;
    }
  }

  // Read in geometry
  const Char_t *detName;
  BrVector3D detOrigin;
  CycDetectorVolume *vol;
  string fn = "DadlGeometry.dat";
  // PathManager *pathMan = PathManager::Instance();
  // const Char_t *dataDir = pathMan->GetDataDir();
  TString fnStr = Form("%s/faustGeom/%s", dataDir, fn.data());
  ifstream cf(fnStr.Data());
  line.ReadLine(cf);
  for (Int_t i = 0; i < 68; i++)
  {
    line.ReadLine(cf);
    Int_t idet;
    Double_t x, y, z, thetaRot, phiRot, psiRot;
    sscanf(line.Data(), "%d %lf %lf %lf %lf %lf %lf", &idet, &x, &y, &z, &thetaRot, &phiRot, &psiRot);
    // printf("idet = %d, x = %f, y = %f, z = %f, theta = %f, phi = %f, psi = %f\n",idet,x,y,z,thetaRot,phiRot,psiRot);

    detOrigin[0] = x;
    detOrigin[1] = y;
    detOrigin[2] = z;
    TGeoRotation rotMatrix(Form("rot%d", idet), phiRot, thetaRot, psiRot);

    detName = Form("Dadl%02d", idet);
    vol = new CycDetectorVolume(detName, detName, 2., 2., 0.000, detOrigin, rotMatrix);
    fDetVol[idet - 1] = vol;
    fFrontPlane[idet - 1] = vol->GetFrontPlane();
  }
  cf.close();

  Clear();
}

void FaustFilter::Begin()
{
  // This method assumes that we are using the run db.
  // Keep in mind that BrDbUpdateModule needs to be in the chain for this
  // to work
  //    CycRunInfoManager *runDb = CycRunInfoManager::Instance();
  //    const CycRunInfo  *run = runDb->GetCurrentRun();
  //    Int_t runNo = run->GetRunNo();
  //    printf("NimrodFilter: Doing begin for run = %d\n",runNo);
  //    Begin(runNo);

  // Set all detectors on for the moment.
}

void FaustFilter::Begin(Int_t run)
{
  fRun = run;
}
void FaustFilter::Clear(Option_t *opt)
{
  for (Int_t i = 0; i < 68; i++)
    fFilteredParticleList[i].Reset();
  fHitList.Clear(); // start on right track
}

void FaustFilter::InitGeometry()
{
  ReadGeometry();
  CalcPlanes();
  // cout << "Initialized FAUST Filter" << endl;
  return;
}

//------------------------------------------------------------------------------

void FaustFilter::ReadGeometry()
{
  // string fn = "FAUST_Corners.txt";
  // string fn = "FAUST_Corners_UnfixedRingE.txt";
  // string fn = "FAUST_Corners.ori";
  // string fn = "FAUSTcornersForClusterRHC.txt";
  string fn = "detnum_corner_z_x_y_r.txt";
  PathManager *pathMan = PathManager::Instance();
  const Char_t *dataDir = pathMan->GetDataDir();
  TString fnStr = Form("%s/faustGeom/%s", dataDir, fn.data());
  // cout << "Reading: " << fnStr.Data() << endl;
  ifstream cf(fnStr.Data());
  string hdr;
  Int_t det, cor;
  Double_t z, x, y, m; //,q;
  if (cf.is_open())
  {
    // cf >> hdr >> hdr >> hdr >> hdr >> hdr >> hdr >> hdr;
    cf >> hdr >> hdr >> hdr >> hdr >> hdr >> hdr >> hdr >> hdr >> hdr >> hdr;
    while (cf.good() && !cf.eof())
    {
      // cf >> det >> cor >> z >> x >> y >> m >> q;
      cf >> det >> cor >> z >> x >> y >> m;
      if (det > -1 && det < 68 && cor >= 1 && cor <= 4)
      {
        fCorner[det][cor - 1][0] = x;
        fCorner[det][cor - 1][1] = y;
        fCorner[det][cor - 1][2] = z;
      }
    }
  }
  cf.close();
  return;
}

//------------------------------------------------------------------------------

void FaustFilter::CalcPlanes()
{
  // cout << "In CalpPlanes" << endl;
  Double_t Aparam, Bparam, Cparam, Dparam;
  for (Int_t det = 0; det < 68; det++)
  {
    // cout << "det = " << det << endl;
    Double_t x1 = fCorner[det][0][0];
    Double_t x2 = fCorner[det][1][0];
    Double_t x3 = fCorner[det][2][0];
    Double_t y1 = fCorner[det][0][1];
    Double_t y2 = fCorner[det][1][1];
    Double_t y3 = fCorner[det][2][1];
    Double_t z1 = fCorner[det][0][2];
    Double_t z2 = fCorner[det][1][2];
    Double_t z3 = fCorner[det][2][2];
    Aparam = (y2 - y1) * (z3 - z1) - (z2 - z1) * (y3 - y1);
    Bparam = (z2 - z1) * (x3 - x1) - (x2 - x1) * (z3 - z1);
    Cparam = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
    Dparam = -x1 * (y2 - y1) * (z3 - z1) + x1 * (z2 - z1) * (y3 - y1) - y1 * (z2 - z1) * (x3 - x1) + y1 * (x2 - x1) * (z3 - z1) - z1 * (x2 - x1) * (y3 - y1) + z1 * (y2 - y1) * (x3 - x1);
    fPlane[det][0] = Aparam;
    fPlane[det][1] = Bparam;
    fPlane[det][2] = Cparam;
    fPlane[det][3] = Dparam;

    /*
        cout << "fCorner 1: " << fCorner[det][0][0] << "\t" << fCorner[det][0][1] << "\t" << fCorner[det][0][2] << endl;
        cout << "fCorner 2: " << fCorner[det][1][0] << "\t" << fCorner[det][1][1] << "\t" << fCorner[det][1][2] << endl;
        cout << "fCorner 3: " << fCorner[det][2][0] << "\t" << fCorner[det][2][1] << "\t" << fCorner[det][2][2] << endl;
        cout << "fCorner 4: " << fCorner[det][3][0] << "\t" << fCorner[det][3][1] << "\t" << fCorner[det][3][2] << endl;
    */

    // Check: how far is the 4th point from the plane defined by the other 3 points?
    Double_t chk_numer = Aparam * fCorner[det][3][0] + Bparam * fCorner[det][3][1] + Cparam * fCorner[det][3][2] + Dparam;
    Double_t chk_denom = sqrt(Aparam * Aparam + Bparam * Bparam + Cparam * Cparam);
    Double_t chk_dist = chk_numer / chk_denom;
    if (fabs(chk_dist) > 5.e-3)
    {
      cout << "det = " << det
           << "\tchk_dist = " << chk_dist
           << endl;
    }

    // Check that all detectors have side length of 2
    // Double_t chk_len01 = sqrt( pow(fCorner[det][0][0]-fCorner[det][1][0],2) + pow(fCorner[det][0][1]-fCorner[det][1][1],2) + pow(fCorner[det][0][2]-fCorner[det][1][2],2) );
    // Double_t chk_len12 = sqrt( pow(fCorner[det][1][0]-fCorner[det][2][0],2) + pow(fCorner[det][1][1]-fCorner[det][2][1],2) + pow(fCorner[det][1][2]-fCorner[det][2][2],2) );
    // Double_t chk_len23 = sqrt( pow(fCorner[det][2][0]-fCorner[det][3][0],2) + pow(fCorner[det][2][1]-fCorner[det][3][1],2) + pow(fCorner[det][2][2]-fCorner[det][3][2],2) );
    // Double_t chk_len30 = sqrt( pow(fCorner[det][3][0]-fCorner[det][0][0],2) + pow(fCorner[det][3][1]-fCorner[det][0][1],2) + pow(fCorner[det][3][2]-fCorner[det][0][2],2) );
    //     if ( fabs(chk_len01-2)>5e-3 ) cout << "det = " << det << "\tchk_len01 = " << chk_len01 << endl;
    //     if ( fabs(chk_len12-2)>5e-3 ) cout << "det = " << det << "\tchk_len12 = " << chk_len12 << endl;
    //     if ( fabs(chk_len23-2)>5e-3 ) cout << "det = " << det << "\tchk_len23 = " << chk_len23 << endl;
    //     if ( fabs(chk_len30-2)>5e-3 ) cout << "det = " << det << "\tchk_len30 = " << chk_len30 << endl;

    // Double_t chk_len02 = sqrt( pow(fCorner[det][0][0]-fCorner[det][2][0],2) + pow(fCorner[det][0][1]-fCorner[det][2][1],2) + pow(fCorner[det][0][2]-fCorner[det][2][2],2) );
    // Double_t chk_len13 = sqrt( pow(fCorner[det][1][0]-fCorner[det][3][0],2) + pow(fCorner[det][1][1]-fCorner[det][3][1],2) + pow(fCorner[det][1][2]-fCorner[det][3][2],2) );
    //     if ( fabs(chk_len02-2*sqrt(2.))>5e-3 ) cout << "det = " << det << "\tchk_len02 = " << chk_len02 << endl;
    //     if ( fabs(chk_len13-2*sqrt(2.))>5e-3 ) cout << "det = " << det << "\tchk_len13 = " << chk_len13 << endl;

  } // end loop over detectors
  return;
}

//==============================================================================

Int_t FaustFilter::GeoFilter(Double_t Q, Double_t F)
{
  if (Q > 90.)
    return -1;
  Double_t interX = 0, interY = 0, interZ = 0;
  // cout << "=====================" << endl;
  // for ( Int_t det=68; det>=1; det-- ) {
  for (Int_t det = 67; det > -1; det--)
  {
    CalcInter(Q, F, det, interX, interY, interZ);
    Double_t area01 = AreaTri(interX, interY, interZ, det, 0, 1);
    Double_t area12 = AreaTri(interX, interY, interZ, det, 1, 2);
    Double_t area23 = AreaTri(interX, interY, interZ, det, 2, 3);
    Double_t area30 = AreaTri(interX, interY, interZ, det, 3, 0);
    Double_t area = area01 + area12 + area23 + area30;

    /*
        cout << "det: " << det << endl;
        cout << "CalcInter: "
             << interX << "\t" << interY << "\t" << interZ
             << endl;
        cout << "a01: " << area01
             << "\ta12: " << area12
             << "\ta23: " << area12
             << "\ta34: " << area12
             << endl;
             cout <<"area: " << area << endl;
    */

    // if ( area <= 4 * 1.0001 ) return det; // the factor 1.0001 allows for rounding errors
    if (area <= 4 * 1.0004)
      return det; // the factor 1.0001 allows for rounding errors
  }
  return -1;
}

Int_t FaustFilter::GeoFilter(Double_t Q, Double_t F, Double_t &interX, Double_t &interY, Double_t &interZ)
{
  if (Q > 90.)
    return -1;
  interX = 0;
  interY = 0;
  interZ = 0;
  // cout << "=====================" << endl;
  // for ( Int_t det=68; det>=1; det-- ) {
  for (Int_t det = 67; det > -1; det--)
  {
    CalcInter(Q, F, det, interX, interY, interZ);

    //     if(det == 60) {
    //        printf("GeoFilter: idet = %d, interX = %f, interY = %f, interZ = %f\n",det,interX,interY,interZ);
    //        }
    Double_t area01 = AreaTri(interX, interY, interZ, det, 0, 1);
    Double_t area12 = AreaTri(interX, interY, interZ, det, 1, 2);
    Double_t area23 = AreaTri(interX, interY, interZ, det, 2, 3);
    Double_t area30 = AreaTri(interX, interY, interZ, det, 3, 0);
    Double_t area = area01 + area12 + area23 + area30;

    /*
        cout << "det: " << det << endl;
        cout << "CalcInter: "
             << interX << "\t" << interY << "\t" << interZ
             << endl;
        cout << "a01: " << area01
             << "\ta12: " << area12
             << "\ta23: " << area12
             << "\ta34: " << area12
             << endl;
    */

    // if ( area <= 4 * 1.0001 ) return det; // the factor 1.0001 allows for rounding errors
    if (area <= 4 * 1.0004)
      return det; // the factor 1.0001 allows for rounding errors
  }
  return -1;
}

//------------------------------------------------------------------------------

void FaustFilter::CalcInter(Double_t Q, Double_t F, Int_t det,
                            Double_t &X, Double_t &Y, Double_t &Z)
{
  Double_t x5 = sin(Q * 3.14159 / 180.) * cos(F * 3.14159 / 180.);
  Double_t y5 = sin(Q * 3.14159 / 180.) * sin(F * 3.14159 / 180.);
  Double_t z5 = cos(Q * 3.14159 / 180.);
  Double_t denom = x5 * fPlane[det][0] + y5 * fPlane[det][1] + z5 * fPlane[det][2];
  Double_t t = -fPlane[det][3] / denom;

/*
  cout <<"x5: " << x5 << endl;
  cout <<"y5: " << y5 << endl;
  cout <<"z5: " << z5 << endl;
  cout <<"fPlane[0]: " << fPlane[det][0] << endl;
  cout <<"fPlane[1]: " << fPlane[det][1] << endl;
  cout <<"fPlane[2]: " << fPlane[det][2] << endl;
  cout <<"fPlane[3]: " << fPlane[det][3] << endl;
  cout <<"t: " << t << endl;
  X = x5 * t;
  */
  Y = y5 * t;
  Z = z5 * t;
  return;
}

//------------------------------------------------------------------------------

Double_t FaustFilter::AreaTri(Double_t x0, Double_t y0, Double_t z0,
                              Int_t det, Int_t c1, Int_t c2)
{
  Double_t x1 = fCorner[det][c1][0];
  Double_t y1 = fCorner[det][c1][1];
  Double_t z1 = fCorner[det][c1][2];
  Double_t x2 = fCorner[det][c2][0];
  Double_t y2 = fCorner[det][c2][1];
  Double_t z2 = fCorner[det][c2][2];
  Double_t dx01 = x1 - x0;
  Double_t dy01 = y1 - y0;
  Double_t dz01 = z1 - z0;
  Double_t dx02 = x2 - x0;
  Double_t dy02 = y2 - y0;
  Double_t dz02 = z2 - z0;
  Double_t dx12 = x2 - x1;
  Double_t dy12 = y2 - y1;
  Double_t dz12 = z2 - z1;
  Double_t a = sqrt(dx01 * dx01 + dy01 * dy01 + dz01 * dz01);
  Double_t b = sqrt(dx02 * dx02 + dy02 * dy02 + dz02 * dz02);
  Double_t c = sqrt(dx12 * dx12 + dy12 * dy12 + dz12 * dz12);
  Double_t s = (a + b + c) / 2;
  Double_t area = sqrt(s * (s - a) * (s - b) * (s - c)); // Heron's formula
  return area;
}

Int_t FaustFilter::Filter(Int_t z, Int_t a, Double_t e, Double_t theta, Double_t phi, TObject *partObj)
{
  if (z < 1)
    return 0; // don't detect neutrons in faust

  //    Double_t interX1, interY1, interZ1;
  //    Int_t detId1 = GeoFilter(theta,phi, interX1, interY1, interZ1);
  //    BrVector3D detectorPos;
  //    if(detId1 > -1) {
  //       CycDetectorVolume *vol = fDetVol[detId1];
  //       BrVector3D intersec;
  //       intersec[0] = interX1;
  //       intersec[1] = interY1;
  //       intersec[2] = interZ1;
  //       detectorPos = vol->GlobalToLocal(intersec,0);
  //       }
  // return detId1;

  Double_t posX = -50000.;
  Double_t posY = -50000.;
  Double_t posZ = -50000.;

  Double_t interX, interY, interZ;

  Int_t detId = IsHit(theta, phi, posX, posY, posZ, interX, interY, interZ);
  if(fIncludeGeometryOnly)
  {
    return detId;
  }

  // if(detId > 0 && detId != detId1+1 && detId1 != -1) printf("detId = %d, detId1 = %d, theta = %f, phi = %f, posX = %f, posX1 = %f, posY = %f, posY1 = %f, posZ = %f, interX = %f, interX1 = %f, interY = %f, interY1 = %f, interZ = %f, interZ1 = %f\n",detId,detId1,theta, phi, posX,detectorPos[0],posY,detectorPos[1],posZ, interX, interX1, interY, interY1, interZ, interZ1);

  if (detId == 0)
    return 0; // no detector hit

  Int_t idet = detId - 1;
  FaustFilterParticle *filteredParticle = &fFilteredParticleList[idet];

  filteredParticle->SetXPos(posX);
  filteredParticle->SetYPos(posY);
  filteredParticle->SetTheta(theta);
  filteredParticle->SetPhi(phi);
  
  filteredParticle->AddHit(idet, partObj);
  // Int_t numHitsInDetector = filteredParticle->GetNumHits();

  Double_t MylarLoss = GetMylarEnergyLoss(idet, e, z, a);
  Double_t E_After_Mylar = e - MylarLoss;

  GetSiEnergyLoss(idet, E_After_Mylar, z, a); // sets global fEnergyLoss variable used in multiple functions so only have to calculate once //add mylar ELoss

  Double_t theta_Smear, phi_Smear;
  // Smear position
  Double_t XPos_Smear = 0;
  Double_t YPos_Smear = 0;
  Double_t ZPos_Smear = 0;
  Double_t interX_Smear = 0;
  Double_t interY_Smear = 0;
  Double_t interZ_Smear = 0;
  Double_t XPos_Rot = 0;
  Double_t YPos_Rot = 0;

  if (fIncludePositionResolution)
  {
   
    GetSmearedPosition(idet, fEnergyLoss, posX, posY, XPos_Smear, YPos_Smear);
    filteredParticle->SetFilteredXPos(XPos_Smear);
    filteredParticle->SetFilteredYPos(YPos_Smear);
     if(fIncludeDetectorPositionSmear)
      {
	//rotate position
	XPos_Rot = XPos_Smear*TMath::Cos(fRotDetSmear[idet]) - YPos_Smear*TMath::Sin(fRotDetSmear[idet]);
	YPos_Rot = YPos_Smear*TMath::Cos(fRotDetSmear[idet]) + XPos_Smear*TMath::Sin(fRotDetSmear[idet]);

	//cout <<"xPos: " << posX << endl;
	//	cout <<"Angle: " << fRotDetSmear[idet] << endl;
	//	cout <<"xRot: " << XPos_Rot << endl;

	//now translate in x and y
	XPos_Smear = XPos_Rot + fXDetSmear[idet];
	YPos_Smear = YPos_Rot + fYDetSmear[idet];
      }
    GetGlobalVector(idet, XPos_Smear, YPos_Smear, ZPos_Smear, interX_Smear, interY_Smear, interZ_Smear);
    GetFaustThetaPhi(interX_Smear, interY_Smear, interZ_Smear, theta_Smear, phi_Smear);
    filteredParticle->SetFilteredTheta(theta_Smear);
    filteredParticle->SetFilteredPhi(phi_Smear);
    
  }
  else
  {
    if(fIncludeDetectorPositionSmear)
      {
	//rotate position
	XPos_Rot = posX*TMath::Cos(fRotDetSmear[idet]) - posY*TMath::Sin(fRotDetSmear[idet]);
	YPos_Rot = posY*TMath::Cos(fRotDetSmear[idet]) + posX*TMath::Sin(fRotDetSmear[idet]);

	//now translate in x and y
	XPos_Smear = XPos_Rot + fXDetSmear[idet];
	YPos_Smear = YPos_Rot + fYDetSmear[idet];
	GetGlobalVector(idet, XPos_Smear, YPos_Smear, ZPos_Smear, interX, interY, interZ);
      }

    GetFaustThetaPhi(interX, interY, interZ, theta, phi);
    filteredParticle->SetFilteredTheta(theta);
    filteredParticle->SetFilteredPhi(phi);
    filteredParticle->SetFilteredXPos(posX);
    filteredParticle->SetFilteredYPos(posY);
  }

  if (fIncludeEnergySmear)
  {
    Double_t SiE_Smear = GetSmearedSiEnergy();
    Double_t CsI_Smear = GetSmearedCsIEnergy(e, MylarLoss);
    Double_t TotalEnergy_Smear = SiE_Smear + CsI_Smear + MylarLoss;
    filteredParticle->SetFilteredE(TotalEnergy_Smear); // smeared 0.8% + 50kev constant for Si and 3% CsI FWHM
    filteredParticle->SetSiEnergy(SiE_Smear);
    filteredParticle->SetCsIEnergy(CsI_Smear);
    // cout <<"SMeared E: " << TotalEnergy_Smear << endl;
  }
  else
  {
    filteredParticle->SetFilteredE(e); // no smearing
    filteredParticle->SetSiEnergy(fEnergyLoss);
    filteredParticle->SetCsIEnergy(e - fEnergyLoss - MylarLoss);
  }

  // if(numHitsInDetector > 1) return 0;  //detector is already hit.

  // coming here means that it cleared the geometry filter.
  Bool_t punchesThrough = PunchesThrough(idet, z, a, e);
  if (punchesThrough) // this means that the particle had enough energy to punch through Si, but not too much to punch out of CsI
  {
    filteredParticle->SetPunchesThrough(kTRUE);
  }
  else
  {
    filteredParticle->SetPunchesThrough(kFALSE);
    if (fIncludePunchthrough)
    {
      return -1;
    }
  }

  Bool_t passThresh = AboveThreshUnderSat(idet, z, a, e, posX, posY);
  if (passThresh)
  {
    filteredParticle->SetAboveThreshUnderSatReject(kFALSE);
  }
  else
  {
    filteredParticle->SetAboveThreshUnderSatReject(kTRUE);
    if (fIncludeThresholdSaturation)
    {
      return -1; // did not pass thresholds of detector
    }
  }
  return detId - 1;
}

Int_t FaustFilter::IsHit(Double_t theta, Double_t phi, Double_t &posX, Double_t &posY, Double_t &posZ, Double_t &interX, Double_t &interY, Double_t &interZ)
{
  const Double_t deg2rad = TMath::Pi() / 180.;

  // convert to radians since that is what routines want!
  theta *= deg2rad;
  phi *= deg2rad;
  Double_t xdir = TMath::Sin(theta) * TMath::Cos(phi);
  Double_t ydir = TMath::Sin(theta) * TMath::Sin(phi);
  Double_t zdir = TMath::Cos(theta);

  BrVector3D origin(0, 0, 0); // target position
  BrVector3D direction(xdir, ydir, zdir);
  BrLine3D line(origin, direction);

  Int_t detId = 0;
  Int_t hitCount = 0;
  // Double_t hitx, hity, hitz;
  // for(Int_t i=0;i<68;i++) {
  for (Int_t i = 68; i > 0; i--)
  {
    Int_t idet = i - 1;
    CycDetectorVolume *vol = fDetVol[idet];
    BrVector3D size = vol->GetSize();
    // Add delta to edges to take into account roundoff error
    Double_t delta = 0.0005;
    Double_t xmin = -size[0] / 2. - delta;
    Double_t xmax = size[0] / 2. + delta;
    Double_t ymin = -size[1] / 2. - delta;
    Double_t ymax = size[1] / 2. + delta;
    BrPlane3D frontPlane = fFrontPlane[idet];

    // Intersection in global coordinates
    BrVector3D intersec = frontPlane.GetIntersectionWithLine(line);

    // Transform to detector local coordinates
    BrVector3D detectorPosition = vol->GlobalToLocal(intersec, 0);
    posX = detectorPosition[0];
    posY = detectorPosition[1];
    posZ = detectorPosition[2]; // this better be 0!
    if (posX > xmin && posX < xmax && posY > ymin && posY < ymax)
    {
      detId = i;
      hitCount++;
      interX = intersec[0];
      interY = intersec[1];
      interZ = intersec[2];
      break;
    }
  }
  // if (hitCount > 1) printf("hitCount = %d\n", hitCount);
  return detId;
}

Bool_t FaustFilter::PunchesThrough(Int_t idet, Int_t z, Int_t a, Double_t energy)
{ // make sure the particle punches through Si to be measured in CsI for PID
  // also make sure it doesnt punch through CsI
  Double_t sithresh = 0.;
  Double_t sicsithresh = 0;
  if (z < 3)
  {
    // for H, He, use precalculated thresholds to optimize time
    Int_t ipar = z - 1 + a - 1;
    sithresh = fSiPunchthrough[idet][ipar];
    sicsithresh = fSiCsIPunchthrough[idet][ipar];
  }
  else
  {
    sithresh = fSrim[idet]->GetPunchthroughEnergy(z, a);
    sicsithresh = 10000000; // set arbitrarily high as Z>2 will never punch through csi
  }
  if (energy < sithresh)
    return kFALSE; // below si punchthrough
  if (energy >= sicsithresh)
    return kFALSE; // above sicsi punchthrough
  else
    return kTRUE; // above si punchthrough and below sicsi punchthrough
}

Bool_t FaustFilter::AboveThreshUnderSat(Int_t idet, Int_t z, Int_t a, Double_t energy, Double_t posX, Double_t posY)
{
  // first grab energy loss through correct thickness of DADL for given isotope
  Double_t SiE = fEnergyLoss;

  // first check X and Y position vs. E saturation. This is a hard cutoff as signal to noise is small
  if (SiE >= TMath::Abs(posX) * fXSatSlope[idet] + fXSatOffset[idet])
  {
    return kFALSE;
  }
  if (SiE >= TMath::Abs(posY) * fYSatSlope[idet] + fYSatOffset[idet])
  {
    return kFALSE;
  }

  // to check low energy threshold, first perform resistive signal splitting with true energy and position
  Double_t F1 = 0;
  Double_t F2 = 0;
  Double_t B1 = 0;
  Double_t B2 = 0;
  SplitSignals(SiE, posX, posY, F1, F2, B1, B2);

  // now add trigger noise to decide if we trigger on the particle (for low DADL E threshold)
  Double_t TrigNoise = gRandom->Gaus(0, fTriggerNoise[idet] / 2.355); // fTriggerNoise is MeV FWHM, argument is sigma
  Double_t F1_T = F1 + TMath::Abs(TrigNoise);                         // adds absolute value of trigger noise, as fast noise can only help you trigger, not hurt you
  Double_t F2_T = F2 + TMath::Abs(TrigNoise);
  Double_t B1_T = B1 + TMath::Abs(TrigNoise);
  Double_t B2_T = B2 + TMath::Abs(TrigNoise);

  // Trigger on each signal
  if (F1_T <= fThreshold[idet])
    return kFALSE;
  if (F2_T <= fThreshold[idet])
    return kFALSE;
  if (B1_T <= fThreshold[idet])
    return kFALSE;
  if (B2_T <= fThreshold[idet])
    return kFALSE;

  else
    return kTRUE;
}

void FaustFilter::SplitSignals(Double_t energyloss, Double_t XPos, Double_t YPos, Double_t &F1, Double_t &F2, Double_t &B1, Double_t &B2)
{

  // first, perform resistive signal splitting with true energy and position
  Double_t FrontOhm = 643; // average effective front face resistance from experimental edge of data on all FAUST DADLs //maybe read in with calibration as this will change exp to exp if changing preamps
  Double_t BackOhm = 322;  // average effective back face resistance from experimental edge of data on all FAUST DADLs //maybe read in with calibration as this will change exp to exp if changing preamps
  Double_t R_Ext = 261;    // actual external resistor value in Ohms

  Double_t R_F1 = ((YPos + 1.0) / 2.0) * FrontOhm + R_Ext;         // resistance experienced for charge going through F1, YPos is scaled to be 0 - 1 instead of -1 to 1
  Double_t R_F2 = (1.0 - ((YPos + 1.0) / 2.0)) * FrontOhm + R_Ext; // resistance experienced for charge going through F2

  Double_t R_B1 = ((XPos + 1.0) / 2.0) * BackOhm + R_Ext;         // resistance experienced for charge going through B1
  Double_t R_B2 = (1.0 - ((XPos + 1.0) / 2.0)) * BackOhm + R_Ext; // resistance experienced for charge going through B2

  Double_t F1_E = R_F2 / (R_F1 + R_F2); // ratio of current through F1 to total
  Double_t B1_E = R_B2 / (R_B1 + R_B2); // ratio of current through B1 to total

  F1 = F1_E * energyloss; // value for F1 in MeV
  F2 = energyloss - F1;
  B1 = B1_E * energyloss;
  B2 = energyloss - B1;

}

void FaustFilter::GetSmearedPosition(Int_t idet, Double_t energyloss, Double_t XPos, Double_t YPos, Double_t &XPos_S, Double_t &YPos_S)
{
  // split signals to get unsmeared f1f2b1b2 values
  Double_t F1 = 0;
  Double_t F2 = 0;
  Double_t B1 = 0;
  Double_t B2 = 0;
  Double_t UncorrNoise = 0;
  SplitSignals(energyloss, XPos, YPos, F1, F2, B1, B2);
  // add noise that contributes to position smearing
  UncorrNoise = gRandom->Gaus(0, fUncorrNoise[idet] / 2.355);
  F1 = F1 + UncorrNoise;
  UncorrNoise = gRandom->Gaus(0, fUncorrNoise[idet] / 2.355);
  F2 = F2 + UncorrNoise;
  UncorrNoise = gRandom->Gaus(0, fUncorrNoise[idet] / 2.355);
  B1 = B1 + UncorrNoise;
  UncorrNoise = gRandom->Gaus(0, fUncorrNoise[idet] / 2.355);
  B2 = B2 + UncorrNoise;

  // calculate position and stretch according to model edge of data using the average resistance values
  Double_t yEdge = 0.55193;
  Double_t xEdge = 0.381517;
  XPos_S = ((B2 - B1) / (B1 + B2)) * (1.0 / xEdge);
  YPos_S = ((F2 - F1) / (F1 + F2)) * (1.0 / yEdge);
}

void FaustFilter::GetSiEnergyLoss(Int_t idet, Double_t energy, Int_t z, Int_t a)
{
  fEnergyLoss = fSrim[idet]->GetEnergyLoss(z, a, energy);
  // cout <<"SiLoss: " << fEnergyLoss << endl;
}

void FaustFilter::GetGlobalVector(Int_t idet, Double_t PosX, Double_t PosY, Double_t PosZ, Double_t &interX, Double_t &interY, Double_t &interZ)
{
  BrVector3D localinter(PosX, PosY, PosZ);
  BrVector3D smearedglobal = fDetVol[idet]->LocalToGlobal(localinter, 0);
  interX = smearedglobal[0];
  interY = smearedglobal[1];
  interZ = smearedglobal[2];
}

Double_t FaustFilter::GetCsIPunchThrough(Int_t idet, Double_t z, Double_t a)
{
  Double_t CsIPunch = fSrimCsI[idet]->GetPunchthroughEnergy(z, a);
  Double_t SiLoss = fSrim[idet]->GetEnergyLossFromResidual(z, a, CsIPunch);
  return CsIPunch + SiLoss;
}

Double_t FaustFilter::GetSmearedSiEnergy()
{
  // emulate constant noise on signals as a constant 50 kev smearing, add to this 0.8% smearing to emulate intrinsic silicon resolution.
  Double_t SiIntrinsic = (fSiSmear/100.) * fEnergyLoss; // % FWHM
  Double_t SiNoise = 0.075;                    // FWHM MeV
  Double_t SiIntrinsicSmear = gRandom->Gaus(0, SiIntrinsic / 2.355);
  Double_t SiNoiseSmear = gRandom->Gaus(0, SiNoise / 2.355);
  return fEnergyLoss + SiIntrinsicSmear + SiNoiseSmear;
}

Double_t FaustFilter::GetSmearedCsIEnergy(Double_t energy, Double_t MylarLoss)
{
  Double_t CsIEnergy = energy - fEnergyLoss - MylarLoss;
  Double_t CsIFWHM = (fCsISmear/100.) * CsIEnergy;
  Double_t CsISmear = gRandom->Gaus(0, CsIFWHM / 2.355);
  if (CsIEnergy <= 1e-3) // do not smear if does not punch through Si
    return 0;

  return CsIEnergy + CsISmear;
}
Double_t FaustFilter::GetMylarEnergyLoss(Int_t idet, Double_t energy, Double_t z, Double_t a)
{
  return fSrimMylar[idet]->GetEnergyLoss(z, a, energy);
}

void FaustFilter::GetFaustThetaPhi(Double_t interX, Double_t interY, Double_t interZ, Double_t &theta, Double_t &phi)
{
  // theta and phi are returned in degrees
  // Assumption that interX, interY, interZ are relative to tgt at 0,0,0

  const Double_t d2rad = TMath::Pi() / 180.;
  const Double_t pi = TMath::Pi();

  Double_t xdir = interX;
  Double_t ydir = interY;
  Double_t zdir = interZ;

  // Add effects of position resolution
  // if (fIncludePositionResolution)
  // {
  //    xdir = gRandom->Gaus(interX, fPositionSigmaX);
  //    ydir = gRandom->Gaus(interY, fPositionSigmaY);
  // }

  Double_t cosTheta = zdir / TMath::Sqrt(xdir * xdir + ydir * ydir + zdir * zdir);

  // Calculate theta and phi from these random directions
  theta = TMath::ACos(cosTheta);

  phi = pi / 2.; // for the case when xdir == 0
  if (ydir < 0.)
    phi = 3. * pi / 2.; // when xdir == 0 and ydir < 0
  if (xdir != 0)
  {
    Double_t tanPhi = ydir / xdir;
    phi = TMath::Abs(TMath::ATan(tanPhi)); // 1st quadrnt
    if (xdir < 0 && ydir > 0)
      phi = pi - phi; // 2nd quadrnt
    else if (xdir < 0 && ydir < 0)
      phi = pi + phi; // 3rd quadrnt
    else if (xdir > 0 && ydir < 0)
      phi = 2. * pi - phi; // 4th quadrnt
  }
  if (ydir == 0.)
    phi = pi;
  if (ydir == 0. && xdir > 0.)
    phi = 2 * pi;

  // Convert to degrees
  theta /= d2rad;
  phi /= d2rad;
}

TObjArray *FaustFilter::GetListOfValidHits()
{
  fHitList.Clear();
  for (Int_t i = 0; i < 68; i++)
  {
    FaustFilterParticle *filteredParticle = &fFilteredParticleList[i];
    if (filteredParticle->ValidHit())
    {
      TObject *partObj = filteredParticle->GetObjectAt(0);
      if (partObj)
        fHitList.Add(partObj);
    }
  }
  return &fHitList;
}

TObjArray *FaustFilter::GetListOfHits()
{
  // return list of all hits

  fHitList.Clear();
  for (Int_t i = 0; i < 68; i++)
  {
    FaustFilterParticle *filteredParticle = &fFilteredParticleList[i];
    Int_t numHits = filteredParticle->GetNumHits();
    if (numHits > 0)
    {
      for (Int_t j = 0; j < numHits; j++)
      {
        TObject *partObj = filteredParticle->GetObjectAt(j);
        if (partObj)
          fHitList.Add(partObj);
      }
    }
  }
  return &fHitList;
}

TObjArray *FaustFilter::GetListOfValidFilteredParticles()
{
  fHitList.Clear();
  for (Int_t i = 0; i < 68; i++)
  {
    FaustFilterParticle *filteredParticle = &fFilteredParticleList[i];
    if (filteredParticle->ValidHit())
      fHitList.Add(filteredParticle);

    // TEST
    // Bool_t validHit = filteredParticle->ValidHit();
    // printf("ListOfValidFilteredParticles:, i = %d, validHit = %d\n",i,validHit);
    // END TEST
  }
  return &fHitList;
}

// $Log: FaustFilter.cxx,v $
// Revision 1.2  2022/11/07 17:45:41  andyhannaman
// Added access functions for energy smearing and to set detector wiggling (x,y,angle) parameters
//
// Revision 1.1  2022/06/20 21:08:31  andyhannaman
// adding Faust filter module to 060821
//
// Revision 1.5  2020/06/03 08:50:41  hagel
// Final modifications to remove warnings thrown by -Wall flag
//
// Revision 1.4  2020/06/03 02:35:01  hagel
// Add -Wall flag
//
// Revision 1.3  2019/09/20 16:10:16  hagel
// Latest
//
// Revision 1.2  2018/09/30 23:13:45  hagel
// Change detId to range from 1 to 68 instead of 0 to 67
//
// Revision 1.1  2018/07/24 00:23:21  hagel
// Initial revision
//
