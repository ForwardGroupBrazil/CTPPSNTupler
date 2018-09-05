void ApplySettings()
{
	lhc_fill = 5279;

	cut1_a = -0.997; cut1_c = -0.781; cut1_si = 0.20;
	cut2_a = -0.972; cut2_c = +0.460; cut2_si = 0.20;
	cut3_a = -1.161; cut3_c = +0.075; cut3_si = 0.15;
	cut4_a = -1.140; cut4_c = +0.025; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(9.6, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(9.0, 15.7);
	selectionRangesX["R_1_N"] = SelectionRange(12.5, 16.0);
	selectionRangesX["R_1_F"] = SelectionRange(12.0, 15.1);

	yAlignmentRangeX["L_1_F"] = SelectionRange(6.0, 12.0);
	yAlignmentRangeX["L_1_N"] = SelectionRange(6.0, 11.0);
	yAlignmentRangeX["R_1_N"] = SelectionRange(5.0, 11.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(4.0, 11.0);
}
