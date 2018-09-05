void ApplySettings()
{
	lhc_fill = 5405;

	cut1_a = -1.002; cut1_c = -0.906; cut1_si = 0.20;
	cut2_a = -1.002; cut2_c = -0.906; cut2_si = 0.20;
	cut3_a = -1.180; cut3_c = -0.637; cut3_si = 0.15;
	cut4_a = -1.180; cut4_c = -0.637; cut4_si = 0.15;

	// TODO
	selectionRangesX["L_1_F"] = SelectionRange(9.8, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(7.1, 14.5);
	selectionRangesX["R_1_N"] = SelectionRange(8.0, 15.5);
	selectionRangesX["R_1_F"] = SelectionRange(7.0, 14.5);
}
