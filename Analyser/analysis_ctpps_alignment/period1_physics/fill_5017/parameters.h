void ApplySettings()
{
	lhc_fill = 5017;

	cut1_a = -1.005; cut1_c = -1.412; cut1_si = 0.20;
	cut2_a = -0.974; cut2_c = +0.472; cut2_si = 0.20;
	cut3_a = -1.160; cut3_c = +0.847; cut3_si = 0.15;
	cut4_a = -1.145; cut4_c = +0.034; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(10.3, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(9.0, 14.5);
	selectionRangesX["R_1_N"] = SelectionRange(8.2, 15.5);
	selectionRangesX["R_1_F"] = SelectionRange(7.8, 14.5);

	yAlignmentRangeX["L_1_F"] = SelectionRange(6.0, 11.0);
	yAlignmentRangeX["L_1_N"] = SelectionRange(6.0, 11.0);
	yAlignmentRangeX["R_1_N"] = SelectionRange(5.0, 11.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(4.0, 11.0);
}
