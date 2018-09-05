void ApplySettings()
{
	alignment_run = true;

	cut1_a = -0.999; cut1_c = +0.046; cut1_si = 0.20;
	cut2_a = -0.963; cut2_c = +0.095; cut2_si = 0.20;
	cut3_a = -1.170; cut3_c = -0.077; cut3_si = 0.15;
	cut4_a = -1.139; cut4_c = +0.400; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(3.5, 11.0);
	selectionRangesX["L_1_N"] = SelectionRange(3.7, 10.5);
	selectionRangesX["R_1_N"] = SelectionRange(2.0, 10.0);
	selectionRangesX["R_1_F"] = SelectionRange(2.0, 10.0);
}
