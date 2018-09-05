void ApplySettings()
{
	lhc_fill = 4976;

	cut1_a = -1.008; cut1_c = -2.285; cut1_si = 0.20;
	cut2_a = -0.971; cut2_c = +0.450; cut2_si = 0.20;
	cut3_a = -1.152; cut3_c = +1.234; cut3_si = 0.15;
	cut4_a = -1.142; cut4_c = +0.110; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(9.8, 16.0);
	selectionRangesX["L_1_N"] = SelectionRange(7.0, 13.5);
	selectionRangesX["R_1_N"] = SelectionRange(8.0, 15.0);
	selectionRangesX["R_1_F"] = SelectionRange(7.0, 14.5);

	yAlignmentRangeX["L_1_F"] = SelectionRange(6.0, 11.0);
	yAlignmentRangeX["L_1_N"] = SelectionRange(6.0, 10.0);
	yAlignmentRangeX["R_1_N"] = SelectionRange(5.0, 11.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(5.0, 11.0);
}
