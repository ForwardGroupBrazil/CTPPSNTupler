void ApplySettings()
{
	lhc_fill = 5261;

	cut1_a = -0.993; cut1_c = -0.820; cut1_si = 0.20;
	cut2_a = -0.970; cut2_c = +0.447; cut2_si = 0.20;
	cut3_a = -1.148; cut3_c = +0.088; cut3_si = 0.15;
	cut4_a = -1.146; cut4_c = +0.030; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(8.7, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(8.3, 15.7);
	selectionRangesX["R_1_N"] = SelectionRange(12.5, 16.0);
	selectionRangesX["R_1_F"] = SelectionRange(12.0, 15.1);

	yAlignmentRangeX["L_1_F"] = SelectionRange(4.0, 11.0);
	yAlignmentRangeX["L_1_N"] = SelectionRange(4.0, 10.0);
	yAlignmentRangeX["R_1_N"] = SelectionRange(4.0, 11.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(4.0, 11.0);
}
