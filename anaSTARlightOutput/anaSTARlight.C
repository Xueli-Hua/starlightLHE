#include "headers.h" 

Double_t shiftDeltaPhi(Double_t dPhi);

const Double_t PI = TMath::Pi();

Int_t oldEvtID;

void anaSTARlight(TString parSpec = "cohJpsi")
{
    TH1::SetDefaultSumw2(kTRUE);

    TH1D *hnEvts        = new TH1D("hnEvts", "hnEvts", 5, -0.5, 4.5);
    TH3D *hMvsPtvsRap_Raw    = new TH3D("hMvsPtvsRap_Raw", "hMvsPtvsRap_Raw; y; p_{T} (GeV); M_{#mu#mu} (GeV); Entries", 50, -2.5, 2.5, 1000, 0, 2, 300, 2, 5);
    TH3D *hMvsPtvsRap_RapSel = new TH3D("hMvsPtvsRap_RapSel", "hMvsPtvsRap_RapSel; y; p_{T} (GeV); M_{#mu#mu} (GeV); Entries", 50, -2.5, 2.5, 1000, 0, 2, 300, 2, 5);
    TH3D *hMvsPtvsRap      = new TH3D("hMvsPtvsRap", "hMvsPtvsRap; y; p_{T} (GeV); M_{#mu#mu} (GeV); Entries", 50, -2.5, 2.5, 1000, 0, 2, 300, 2, 5);
    TH1D *hPt2_RapSel      = new TH1D("hPt2_RapSel", "hPt2_RapSel; p_{T}^{2} (GeV)^{2}", 500, 0, 0.1);
    TH3D *hPosPhivsEtavsPt = new TH3D("hPosPhivsEtavsPt", "hPosPhivsEtavsPt; p_{T} (GeV); #eta; #phi; Entries", 500, 0, 5, 100, -2.5, 2.5, 120, -PI, PI);
    TH3D *hNegPhivsEtavsPt = new TH3D("hNegPhivsEtavsPt", "hNegPhivsEtavsPt; p_{T} (GeV); #eta; #phi; Entries", 500, 0, 5, 100, -2.5, 2.5, 120, -PI, PI);
    TH2D *hNegPtvsPosPt    = new TH2D("hNegPtvsPosPt", "hNegPtvsPosPt; #mu^{+} p_{T} (GeV); #mu^{-} p_{T} (GeV)", 500, 0, 5, 500, 0, 5);
    TH2D *hNegEtavsPosEta  = new TH2D("hNegEtavsPosEta", "hNegEtavsPosEta; #mu^{+} #eta; #mu^{-} #eta", 100, -2.5, 2.5, 100, -2.5, 2.5);
    TH2D *hNegPhivsPosPhi  = new TH2D("hNegPhivsPosPhi", "hNegPhivsPosPhi; #mu^{+} #phi; #mu^{-} #phi", 120, -PI, PI, 120, -PI, PI);
    TH3D *hPosPhivsEtavsP = new TH3D("hPosPhivsEtavsP", "hPosPhivsEtavsP; p (GeV); #eta; #phi; Entries", 500, 0, 5, 100, -2.5, 2.5, 120, -PI, PI);
    TH3D *hNegPhivsEtavsP = new TH3D("hNegPhivsEtavsP", "hNegPhivsEtavsP; p (GeV); #eta; #phi; Entries", 500, 0, 5, 100, -2.5, 2.5, 120, -PI, PI);

    ifstream infile(Form("%s.out", parSpec.Data()));

    string temp_string, temp;
    istringstream curstring;

    int nEvts = 0;

    TLorentzVector posParMom, negParMom;
    Int_t nPosMu=0, nNegMu=0;
    while (getline(infile,temp_string)) {
        curstring.clear(); // needed when using several tims istringstream::str(string)
        curstring.str(temp_string);

        if(strstr(temp_string.c_str(),"EVENT")){
            nEvts++;
            posParMom.SetPtEtaPhiE(0, 0, 0, 0);
            negParMom.SetPtEtaPhiE(0, 0, 0, 0);
            nPosMu = 0;
            nNegMu = 0;
        }
        else if(strstr(temp_string.c_str(),"TRACK")){
            int useless, part_n, geant3_id, pdg_id;
            double px, py, pz;
            curstring >> temp >> geant3_id >> px >> py >> pz >> part_n >> useless >> useless >> pdg_id;
            //cout<< px << "   " << py << "   " << pz << "   " << pdg_id << endl;

            if(TMath::Abs(pdg_id) != 13) continue;

            TParticle par (pdg_id, 0, 0, 0, 0, 0, px, py, pz, 0.0, 0.0, 0.0, 0.0, 0.0);
            double mass = par.GetMass();

            float E = TMath::Sqrt(pow(mass,2)+pow(par.P(),2));

            if(pdg_id == 13){
                negParMom.SetPtEtaPhiE(par.Pt(), par.Eta(), par.Phi(), E);
                nNegMu++;
            }
            else if(pdg_id == -13 ){
                posParMom.SetPtEtaPhiE(par.Pt(), par.Eta(), par.Phi(), E);
                nPosMu++;
            }
            else  cout<<"particle information is not correct"<<endl;
        }

        if(nPosMu==1 && nNegMu==1){
            if(nEvts%100000 == 0) cout<<"Working on "<< nEvts/100000 <<"-th 100k event ..."<<endl;

            TLorentzVector motherPar = negParMom + posParMom;
            double y = motherPar.Rapidity();

            hnEvts->Fill(0);
            hMvsPtvsRap_Raw->Fill(motherPar.Rapidity(), motherPar.Pt(), motherPar.M());

            if(TMath::Abs(y)>2.4) continue;

            hnEvts->Fill(1);
            hMvsPtvsRap_RapSel->Fill(motherPar.Rapidity(), motherPar.Pt(), motherPar.M());
            hPt2_RapSel->Fill(motherPar.Pt()*motherPar.Pt());

            if(
                    TMath::Abs(posParMom.Eta())>1.3 && TMath::Abs(posParMom.Eta())<2.4
                    && TMath::Abs(negParMom.Eta())>1.3 && TMath::Abs(negParMom.Eta())<2.4
              ){
                hnEvts->Fill(2);

                hMvsPtvsRap->Fill(motherPar.Rapidity(), motherPar.Pt(), motherPar.M());

                hPosPhivsEtavsPt->Fill(posParMom.Pt(), posParMom.Eta(), posParMom.Phi());
                hNegPhivsEtavsPt->Fill(negParMom.Pt(), negParMom.Eta(), negParMom.Phi());

                hNegPtvsPosPt->Fill(posParMom.Pt(), negParMom.Pt());
                hNegEtavsPosEta->Fill(posParMom.Eta(), negParMom.Eta());
                hNegPhivsPosPhi->Fill(posParMom.Phi(), negParMom.Phi());
            }
            if(motherPar.Pt()>0.2) continue;
            hPosPhivsEtavsP->Fill(posParMom.P(), posParMom.Eta(), posParMom.Phi());
            hNegPhivsEtavsP->Fill(negParMom.P(), negParMom.Eta(), negParMom.Phi());
        }

    } // reading loop of the input file

    TFile *fOut = new TFile(Form("%s.root", parSpec.Data()), "recreate");
    fOut->cd();
    hnEvts->Write();
    hMvsPtvsRap_Raw->Write();
    hMvsPtvsRap_RapSel->Write();
    hPt2_RapSel->Write();
    hMvsPtvsRap->Write();
    hPosPhivsEtavsPt->Write();
    hNegPhivsEtavsPt->Write();
    hNegPtvsPosPt->Write();
    hNegEtavsPosEta->Write();
    hNegPhivsPosPhi->Write();

    hPosPhivsEtavsP->Write();
    hNegPhivsEtavsP->Write();
    fOut->Close();
}

Double_t shiftDeltaPhi(Double_t dPhi)
{
    if (isnan(dPhi))
        return -999;

    while (dPhi < -PI)
        dPhi += 2 * PI;
    while (dPhi >= PI)
        dPhi -= 2 * PI;

    return dPhi;
}
