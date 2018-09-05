all: DarkMatterSearchLeadingMuons DarkMatterSearchLeadingMuonsJets DarkMatterSearchLeadingElectrons DarkMatterSearchLeadingElectronsJets DarkMatterSearchLeadingMuonsFinn

DarkMatterSearchLeadingMuons: DarkMatterSearchLeadingMuons.cc
	g++ --std=c++11 `root-config --libs` `root-config --cflags` \
        DarkMatterSearchLeadingMuons.cc -o DarkMatterSearchLeadingMuons

DarkMatterSearchLeadingMuonsJets: DarkMatterSearchLeadingMuonsJets.cc
	g++ --std=c++11 `root-config --libs` `root-config --cflags` \
        DarkMatterSearchLeadingMuonsJets.cc -o DarkMatterSearchLeadingMuonsJets

DarkMatterSearchLeadingElectrons: DarkMatterSearchLeadingElectrons.cc
	g++ --std=c++11 `root-config --libs` `root-config --cflags` \
        DarkMatterSearchLeadingElectrons.cc -o DarkMatterSearchLeadingElectrons

DarkMatterSearchLeadingElectronsJets: DarkMatterSearchLeadingElectronsJets.cc
	g++ --std=c++11 `root-config --libs` `root-config --cflags` \
        DarkMatterSearchLeadingElectronsJets.cc -o DarkMatterSearchLeadingElectronsJets

DarkMatterSearchLeadingMuonsFinn: DarkMatterSearchLeadingMuonsFinn.cc
	g++ --std=c++11 `root-config --libs` `root-config --cflags` \
	DarkMatterSearchLeadingMuonsFinn.cc -o DarkMatterSearchLeadingMuonsFinn

clean:
	rm DarkMatterSearchLeadingMuons
	rm DarkMatterSearchLeadingMuonsJets
	rm DarkMatterSearchLeadingElectrons
	rm DarkMatterSearchLeadingElectronsJets
	rm DarkMatterSearchLeadingMuonsFinn


