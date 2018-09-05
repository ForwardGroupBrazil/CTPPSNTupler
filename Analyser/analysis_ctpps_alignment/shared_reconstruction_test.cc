#include "common.h"

#include "shared_track.h"
#include "shared_alignment.h"
#include "shared_reconstruction.h"
#include "shared_fill_info.h"

#include "input_files.h"
#include "parameters.h"

#include "stat.h"

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

struct PlotGroup
{
	map<unsigned int, TH1D*> h_x;

	map<unsigned int, TProfile*> p_y_vs_x;

	map<unsigned int, TH1D*> h_th_x;
	map<unsigned int, TH1D*> h_th_y;

	map<unsigned int, TH1D*> h_xi;


	PlotGroup()
	{
		for (int rpId : {3, 2, 102, 103})
		{
			h_x[rpId] = new TH1D("", ";x", 200, 0., 20.);

			p_y_vs_x[rpId] = new TProfile("", ";x;y", 200, 0., 20.);

			h_th_x[rpId] = new TH1D("", ";th_x", 200, -0.02, +0.02);
			h_th_y[rpId] = new TH1D("", ";th_y", 200, -0.02, +0.02);

			h_xi[rpId] = new TH1D("", ";xi", 200, 0., 0.2);
		}

	}

	void Write() const
	{
		char buf[100];

		for (const auto &it : h_x)
		{
			sprintf(buf, "h_x_%u", it.first);
			it.second->Write(buf);
		}

		for (const auto &it : p_y_vs_x)
		{
			sprintf(buf, "p_y_vs_x_%u", it.first);
			it.second->Write(buf);
		}

		for (const auto &it : h_th_x)
		{
			sprintf(buf, "h_th_x_%u", it.first);
			it.second->Write(buf);
		}

		for (const auto &it : h_th_y)
		{
			sprintf(buf, "h_th_y_%u", it.first);
			it.second->Write(buf);
		}

		for (const auto &it : h_xi)
		{
			sprintf(buf, "h_xi_%u", it.first);
			it.second->Write(buf);
		}
	}
};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int main()
{
	// apply settings
	ApplySettings();

	InitReconstruction();
	InitFillInfoCollection();

	// get input
	InitInputFiles();
	fwlite::ChainEvent event(input_files);

	printf("* events in input chain: %llu\n", event.size());

	// ouput file
	TFile *f_out = new TFile("shared_reconstruction_test.root", "recreate");

	// book histograms
	// map: method idx, cuts on(1)/off(0) --> plots
	map<unsigned int, PlotGroup> plots;

	//  load alignment collection
	AlignmentResultsCollection alignmentCollection;
	alignmentCollection.Load("../../shared_alignment/collect_alignments.out");

	// loop over the chain entries
	unsigned long int ev_count = 0;
	unsigned long int tr_sel_count = 0;

	AlignmentResults *alignments = NULL;
	unsigned int prev_run = 0;

	for (event.toBegin(); ! event.atEnd(); ++event)
	{
		// TODO: comment out
		/*
		if (ev_count > 1000)
			break;
		*/

		ev_count++;

		// get track data
		fwlite::Handle< DetSetVector<TotemRPLocalTrack> > tracks;
		tracks.getByLabel(event, "totemRPLocalTrackFitter");

		TrackDataCollection tr;

		for (const auto &ds : *tracks)
		{
			const auto &rpId = ds.detId();

			for (const auto &track : ds)
			{
				if (rpId == 2 || rpId == 3)
					tr[rpId] = track;

				if (rpId == 102 || rpId == 103)
					tr[rpId] = track;
			}
		}

		// cut evaluation
		double cq1 = tr[3].x + cut1_a * tr[2].x + cut1_c;
		bool cut1_val = (fabs(cq1) < n_si * cut1_si);

		double cq2 = tr[103].x + cut2_a * tr[102].x + cut2_c;
		bool cut2_val = (fabs(cq2) < n_si * cut2_si);

		double cq3 = tr[3].y + cut3_a * tr[2].y + cut3_c;
		bool cut3_val = (fabs(cq3) < n_si * cut3_si);

		double cq4 = tr[103].y + cut4_a * tr[102].y + cut4_c;
		bool cut4_val = (fabs(cq4) < n_si * cut4_si);
		
		bool all_tracks_L = tr[2].valid && tr[3].valid;
		bool all_tracks_R = tr[102].valid && tr[103].valid;

		bool cuts_L = all_tracks_L && (!cut1_apply || cut1_val) && (!cut3_apply || cut3_val);
		bool cuts_R = all_tracks_R && (!cut2_apply || cut2_val) && (!cut4_apply || cut4_val);

		// increase counter
		if (cuts_L || cuts_R)
			tr_sel_count++;

		// get alignment
		if (prev_run != event.id().run() || alignments == NULL)
		{
			const auto &fillInfo = fillInfoCollection.FindByRun(event.id().run());
			prev_run = event.id().run();

			const auto alignment_it = alignmentCollection.find(fillInfo.alignmentTag);
			if (alignment_it == alignmentCollection.end())
			{
				printf("ERROR: no alignment for tag '%s'.\n", fillInfo.alignmentTag.c_str());
				return 1;
			}

			printf("INFO: loaded alignment with tag '%s'.\n", fillInfo.alignmentTag.c_str());

			alignments = &alignment_it->second;
		}

		// apply alignments
		TrackDataCollection tr_al = alignments->Apply(tr);

		// loop over RPs
		map<unsigned int, ProtonData> perRPProtonReco;

		for (unsigned int rpId : {2, 3, 102, 103})
		{
			bool left = ((rpId / 100) == 0);
			bool cuts = (left) ? cuts_L : cuts_R;

			const auto &track = tr_al[rpId];
			if (!track.valid)
				continue;

			ProtonData proton = ReconstructProton({{rpId, track}}, left);

			perRPProtonReco[rpId] = proton;
			
			plots[0].h_x[rpId]->Fill(track.x);
			plots[0].p_y_vs_x[rpId]->Fill(track.x, track.y);
			plots[0].h_th_x[rpId]->Fill(track.th_x);
			plots[0].h_th_y[rpId]->Fill(track.th_y);

			if (cuts)
			{
				plots[1].h_x[rpId]->Fill(track.x);
				plots[1].p_y_vs_x[rpId]->Fill(track.x, track.y);
				plots[1].h_th_x[rpId]->Fill(track.th_x);
				plots[1].h_th_y[rpId]->Fill(track.th_y);
			}

			if (proton.valid)
			{
				plots[0].h_xi[rpId]->Fill(proton.xi);

				if (cuts)
					plots[1].h_xi[rpId]->Fill(proton.xi);
			}
		}
	}

	printf("* events processed: %lu\n", ev_count);
	printf("* tracks selected: %lu\n", tr_sel_count);

	// save plots
	for (const auto &cit : plots)
	{
		TDirectory *c_dir = f_out->mkdir( (cit.first) ? "with cuts" : "without cuts" );

		gDirectory = c_dir;

		cit.second.Write();
	}

	// clean up
	delete f_out;
	return 0;
}
