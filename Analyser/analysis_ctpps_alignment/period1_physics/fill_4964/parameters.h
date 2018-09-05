void ApplySettings()
{
	lhc_fill = 4964;

	cut1_a = -1.005; cut1_c = -2.357; cut1_si = 0.20;
	cut2_a = -0.970; cut2_c = +0.434; cut2_si = 0.20;
	cut3_a = -1.160; cut3_c = +1.155; cut3_si = 0.15;
	cut4_a = -1.145; cut4_c = +0.034; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(9.8, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(7.1, 14.5);
	selectionRangesX["R_1_N"] = SelectionRange(8.0, 15.5);
	selectionRangesX["R_1_F"] = SelectionRange(7.0, 14.5);

	yAlignmentRangeX["L_1_F"] = SelectionRange(6.0, 11.0);
	yAlignmentRangeX["L_1_N"] = SelectionRange(5.5, 10.5);
	yAlignmentRangeX["R_1_N"] = SelectionRange(5.0, 10.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(4.0, 10.0);
}
