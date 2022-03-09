#include <TString.h>
#include <TROOT.h>

void LoadVtx(TString dname)
{
gROOT->Reset();
TString Action = ".L SolGeom" + dname + ".cc+";
gROOT->ProcessLine(Action);
gROOT->ProcessLine(".L TrkUtil.cc+");
gROOT->ProcessLine(".L SolTrack.cc+");
}
