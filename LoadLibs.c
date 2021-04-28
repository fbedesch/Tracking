#include <TString.h>
#include <TROOT.h>

void LoadLibs(TString dname)
{
gROOT->Reset();
TString Action = ".L SolGeom" + dname + ".cc+";
gROOT->ProcessLine(Action);
gROOT->ProcessLine(".L SolTrack.cc+");
gROOT->ProcessLine(".L AcceptanceClx.cc+");
gROOT->ProcessLine(".L SolGridCov.cc+");
gROOT->ProcessLine(".L TrkUtil.cc+");
gROOT->ProcessLine(".L ObsTrk.cc+");
gROOT->ProcessLine(".L VertexFit.cc+");
gROOT->ProcessLine(".include $PYTHIA8/include");
gROOT->ProcessLine(".L $PYTHIA8/lib/libpythia8.so");
}
