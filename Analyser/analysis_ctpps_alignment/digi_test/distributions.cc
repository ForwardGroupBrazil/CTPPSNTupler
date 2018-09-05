#include "input_files.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"
#include "DataFormats/CTPPSReco/interface/TotemRPUVPattern.h"
#include "DataFormats/CTPPSReco/interface/TotemRPCluster.h"
#include "DataFormats/TotemDigi/interface/TotemRPDigi.h"
#include "DataFormats/TotemDigi/interface/TotemVFATStatus.h"
#include "DataFormats/TotemRPDetId/interface/TotemRPDetId.h"

#include <vector>
#include <string>

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

const unsigned short no_of_strips_ = 512; 
const double pitch_ = 66E-3;
const double y_width_ = 36.07;
const double last_strip_to_border_dist_ = 1.4175;

//----------------------------------------------------------------------------------------------------

struct SensorGeometry
{
	double c_x;		// mm
	double c_y;		// mm
	double c_z;		// m
	double rod_x;
	double rod_y;
};

map<unsigned int, SensorGeometry> sensorGeometry;


struct RPGeometry
{
	double c_x;	// mm
	double c_y;	// mm
	double c_z;	// m
	double rod_U_x;
	double rod_U_y;
	double rod_V_x;
	double rod_V_y;
};

map<unsigned int, RPGeometry> rpGeometry;

//----------------------------------------------------------------------------------------------------

struct RPPlots
{
	unsigned int id;
	
	vector<TH1D *> v_plane_h_digi;		// in strip number
	vector<TH1D *> v_plane_h_digi_abs;	// in mm

	TH2D *h2_digi_corr_x_y;
	TH2D *h2_digi_w_corr_x_y;

	TH2D *h2_tracks_x_y = NULL;
	TH1D *h_tracks_U = NULL;
	TH1D *h_tracks_V = NULL;

	RPPlots(unsigned int _id) : id(_id)
	{
		for (unsigned int i = 0; i < 10; ++i)
		{
			v_plane_h_digi.push_back(new TH1D("", "strip", 512, -0.5, 511.5));
			v_plane_h_digi_abs.push_back(new TH1D("", "mm", 500, -100. * 66E-3, +400. * 66E-3));
		}

		h2_digi_corr_x_y = new TH2D("", ";x;y", 100, 0., 20., 100, -10., +10.);
		h2_digi_w_corr_x_y = new TH2D("", ";x;y", 100, 0., 20., 100, -10., +10.);

		h2_tracks_x_y = new TH2D("", ";x;y", 100, 0., 20., 100, -10., +10.);

		h_tracks_U = new TH1D("", "", 500, -100. * 66E-3, +400. * 66E-3);
		h_tracks_V = new TH1D("", "", 500, -100. * 66E-3, +400. * 66E-3);
	}

	void Process(const DetSetVector<TotemRPDigi> &digis, const DetSetVector<TotemRPLocalTrack> &tracks)
	{
		unsigned int N_pl0 = 0, N_pl1 = 0;

		for (const auto &drp : digis)
		{
			unsigned int decId = TotemRPDetId::rawToDecId(drp.detId());
			unsigned int rpId = decId / 10;
			unsigned int plIdx = decId % 10;

			if (rpId != id)
				continue;

			const auto geom = sensorGeometry[decId];

			for (const auto &d : drp)
			{
				v_plane_h_digi[plIdx]->Fill(d.getStripNumber());

				double strip_no = d.getStripNumber();

				double reco = last_strip_to_border_dist_ + (no_of_strips_-1)*pitch_ - y_width_/2. - strip_no * pitch_;
				double offset = geom.c_x * geom.rod_x + geom.c_y * geom.rod_y;

				v_plane_h_digi_abs[plIdx]->Fill(reco + offset);

				if (plIdx == 0)
					N_pl0++;

				if (plIdx == 1)
					N_pl1++;
			}
		}

		unsigned int plane_u = 10*id + 1;
		unsigned int plane_v = 10*id + 0;

		const auto &geom_u = sensorGeometry[plane_u];
		const auto &geom_v = sensorGeometry[plane_v];
		
		double offset_u = geom_u.c_x * geom_u.rod_x + geom_u.c_y * geom_u.rod_y;
		double offset_v = geom_v.c_x * geom_v.rod_x + geom_v.c_y * geom_v.rod_y;

		double det = geom_u.rod_x * geom_v.rod_y - geom_u.rod_y * geom_v.rod_x;

		for (const auto &dse_u : digis)
		{
			unsigned int decId_u = TotemRPDetId::rawToDecId(dse_u.detId());
			if (decId_u != plane_u)
				continue;

			for (const auto &dse_v : digis)
			{
				unsigned int decId_v = TotemRPDetId::rawToDecId(dse_v.detId());
				if (decId_v != plane_v)
					continue;

				for (const auto &d_u : dse_u)
				{
					for (const auto &d_v : dse_v)
					{
						double reco_u = last_strip_to_border_dist_ + (no_of_strips_-1)*pitch_ - y_width_/2. - d_u.getStripNumber() * pitch_;
						double reco_v = last_strip_to_border_dist_ + (no_of_strips_-1)*pitch_ - y_width_/2. - d_v.getStripNumber() * pitch_;

						double u = reco_u + offset_u;
						double v = reco_v + offset_v;

						double x = (geom_v.rod_y * u - geom_u.rod_y * v) / det;
						double y = (- geom_v.rod_x * u + geom_u.rod_x * v) / det;

						double w = 1. / N_pl0 / N_pl1;

						h2_digi_corr_x_y->Fill(x, y);
						h2_digi_w_corr_x_y->Fill(x, y, w);
					}
				}
			}
		}

		for (const auto &trp : tracks)
		{
			if (trp.detId() != id)
				continue;

			const auto geom = rpGeometry[id];

			for (const auto &tr : trp)
			{
				h2_tracks_x_y->Fill(tr.getX0(), tr.getY0());

				double U = tr.getX0() * geom.rod_U_x + tr.getY0() * geom.rod_U_y;
				double V = tr.getX0() * geom.rod_V_x + tr.getY0() * geom.rod_V_y;

				h_tracks_U->Fill(U);
				h_tracks_V->Fill(V);
			}
		}
	}

	void Write() const
	{
		TDirectory *topDir = gDirectory;

		char buf[100];
		sprintf(buf, "RP %u", id);
		gDirectory = topDir->mkdir(buf);

		for (unsigned int i = 0; i < 10; ++i)
		{
			sprintf(buf, "h_plane_%u_digi", i);
			v_plane_h_digi[i]->Write(buf);

			sprintf(buf, "h_plane_%u_digi_abs", i);
			v_plane_h_digi_abs[i]->Write(buf);
		}
		
		h2_digi_corr_x_y->Write("h2_digi_corr_x_y");
		h2_digi_w_corr_x_y->Write("h2_digi_w_corr_x_y");

		h2_tracks_x_y->Write("h2_tracks_x_y");
		
		h_tracks_U->Write("h_tracks_U");
		h_tracks_V->Write("h_tracks_V");

		gDirectory = topDir;
	}
};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int main()
{
	// sensor geometry
	sensorGeometry[  20] = {    17.986  ,    -0.457  ,  -203.8067 ,    0.702 ,   -0.712 };
	sensorGeometry[  21] = {    17.993  ,    -1.400  ,  -203.8113 ,    0.712 ,    0.702 };
	sensorGeometry[  22] = {    17.969  ,    -0.441  ,  -203.8157 ,    0.702 ,   -0.712 };
	sensorGeometry[  23] = {    17.976  ,    -1.418  ,  -203.8203 ,    0.712 ,    0.702 };
	sensorGeometry[  24] = {    17.960  ,    -0.432  ,  -203.8247 ,    0.701 ,   -0.713 };
	sensorGeometry[  25] = {    17.939  ,    -1.457  ,  -203.8293 ,    0.712 ,    0.702 };
	sensorGeometry[  26] = {    17.930  ,    -0.403  ,  -203.8337 ,    0.702 ,   -0.713 };
	sensorGeometry[  27] = {    17.914  ,    -1.482  ,  -203.8383 ,    0.713 ,    0.701 };
	sensorGeometry[  28] = {    17.893  ,    -0.366  ,  -203.8427 ,    0.701 ,   -0.713 };
	sensorGeometry[  29] = {    17.916  ,    -1.480  ,  -203.8473 ,    0.713 ,    0.701 };

	sensorGeometry[  30] = {    15.430  ,     3.311  ,  -212.5297 ,    0.782 ,   -0.624 };
	sensorGeometry[  31] = {    15.597  ,     2.239  ,  -212.5343 ,    0.624 ,    0.781 };
	sensorGeometry[  32] = {    15.439  ,     3.305  ,  -212.5387 ,    0.781 ,   -0.624 };
	sensorGeometry[  33] = {    15.608  ,     2.253  ,  -212.5433 ,    0.624 ,    0.781 };
	sensorGeometry[  34] = {    15.488  ,     3.267  ,  -212.5477 ,    0.781 ,   -0.624 };
	sensorGeometry[  35] = {    15.606  ,     2.251  ,  -212.5523 ,    0.624 ,    0.781 };
	sensorGeometry[  36] = {    15.474  ,     3.278  ,  -212.5567 ,    0.781 ,   -0.624 };
	sensorGeometry[  37] = {    15.621  ,     2.270  ,  -212.5613 ,    0.624 ,    0.781 };
	sensorGeometry[  38] = {    15.528  ,     3.237  ,  -212.5657 ,    0.781 ,   -0.625 };
	sensorGeometry[  39] = {    15.642  ,     2.297  ,  -212.5703 ,    0.625 ,    0.781 };

	sensorGeometry[1020] = {    16.121  ,     0.220  ,   203.8067 ,    0.695 ,    0.719 };
	sensorGeometry[1021] = {    16.191  ,     1.208  ,   203.8113 ,    0.720 ,   -0.694 };
	sensorGeometry[1022] = {    16.139  ,     0.238  ,   203.8157 ,    0.694 ,    0.720 };
	sensorGeometry[1023] = {    16.159  ,     1.239  ,   203.8203 ,    0.720 ,   -0.694 };
	sensorGeometry[1024] = {    16.132  ,     0.231  ,   203.8247 ,    0.695 ,    0.719 };
	sensorGeometry[1025] = {    16.132  ,     1.266  ,   203.8293 ,    0.720 ,   -0.694 };
	sensorGeometry[1026] = {    16.134  ,     0.233  ,   203.8337 ,    0.694 ,    0.720 };
	sensorGeometry[1027] = {    16.110  ,     1.288  ,   203.8383 ,    0.720 ,   -0.694 };
	sensorGeometry[1028] = {    16.153  ,     0.253  ,   203.8427 ,    0.694 ,    0.720 };
	sensorGeometry[1029] = {    16.086  ,     1.312  ,   203.8473 ,    0.720 ,   -0.694 };

	sensorGeometry[1030] = {    15.601  ,     2.502  ,   212.5297 ,    0.586 ,    0.810 };
	sensorGeometry[1031] = {    15.600  ,     3.522  ,   212.5343 ,    0.810 ,   -0.586 };
	sensorGeometry[1032] = {    15.636  ,     2.548  ,   212.5387 ,    0.586 ,    0.810 };
	sensorGeometry[1033] = {    15.557  ,     3.554  ,   212.5433 ,    0.810 ,   -0.587 };
	sensorGeometry[1034] = {    15.653  ,     2.571  ,   212.5477 ,    0.586 ,    0.810 };
	sensorGeometry[1035] = {    15.512  ,     3.589  ,   212.5523 ,    0.810 ,   -0.586 };
	sensorGeometry[1036] = {    15.668  ,     2.591  ,   212.5567 ,    0.585 ,    0.811 };
	sensorGeometry[1037] = {    15.461  ,     3.628  ,   212.5613 ,    0.810 ,   -0.586 };
	sensorGeometry[1038] = {    15.711  ,     2.648  ,   212.5657 ,    0.585 ,    0.811 };
	sensorGeometry[1039] = {    15.424  ,     3.655  ,   212.5703 ,    0.811 ,   -0.585 };

	// RP geometry
	rpGeometry[  2] = {    17.948  ,    -0.934  ,  -203.8270 ,    0.713 ,    0.702 ,    0.702 ,   -0.713 };
	rpGeometry[  3] = {    15.543  ,     2.771  ,  -212.5500 ,    0.624 ,    0.781 ,    0.781 ,   -0.624 };
	rpGeometry[102] = {    16.136  ,     0.749  ,   203.8270 ,    0.720 ,   -0.694 ,    0.695 ,    0.719 };
	rpGeometry[103] = {    15.582  ,     3.081  ,   212.5500 ,    0.810 ,   -0.586 ,    0.586 ,    0.810 };

	// get input
	InitInputFiles();

	//input_files.clear();
	//input_files.push_back("reco_test.root");
	//input_files.push_back("../test3/reco.root");

	if (input_files.empty())
	{
		printf("WARNING: input_files empty.\n");
		return 0;
	}

	for (const auto &f : input_files)
		printf("+ %s\n", f.c_str());

	fwlite::ChainEvent ev(input_files);

	printf("* events in input chain: %llu\n", ev.size());

	// ouput file
	TFile *f_out = new TFile("distributions.root", "recreate");

	// book histograms
	vector<RPPlots> rpPlots;
	rpPlots.push_back(RPPlots(3));
	rpPlots.push_back(RPPlots(2));
	rpPlots.push_back(RPPlots(102));
	rpPlots.push_back(RPPlots(103));

	// loop over the chain entries
	unsigned long int ev_count = 0;
	for (ev.toBegin(); ! ev.atEnd(); ++ev)
	{
		// TODO: comment out
		if (ev_count > 100000)
			break;

		ev_count++;

		// get event data
		fwlite::Handle< DetSetVector<TotemRPDigi> > digis;
		digis.getByLabel(ev, "totemRPRawToDigi", "RP");
		
		fwlite::Handle< DetSetVector<TotemRPLocalTrack> > tracks;
		tracks.getByLabel(ev, "totemRPLocalTrackFitter");

		for (auto &p : rpPlots)
			p.Process(*digis, *tracks);
	}

	printf("* events processed: %lu\n", ev_count);

	// save histograms
	gDirectory = f_out;

	// save histograms
	for (const auto &p : rpPlots)
		p.Write();

	// clean up
	delete f_out;
	return 0;
}
