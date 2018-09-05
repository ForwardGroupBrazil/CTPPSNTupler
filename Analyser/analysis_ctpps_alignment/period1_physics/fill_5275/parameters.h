void ApplySettings()
{
	lhc_fill = 5275;

	cut1_a = -0.995; cut1_c = -0.785; cut1_si = 0.20;
	cut2_a = -0.970; cut2_c = +0.429; cut2_si = 0.20;
	cut3_a = -1.157; cut3_c = +0.093; cut3_si = 0.15;
	cut4_a = -1.138; cut4_c = +0.029; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(9.5, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(8.5, 15.7);
	selectionRangesX["R_1_N"] = SelectionRange(12.5, 16.0);
	selectionRangesX["R_1_F"] = SelectionRange(12.0, 15.1);

	yAlignmentRangeX["L_1_F"] = SelectionRange(5.0, 11.0);
	yAlignmentRangeX["L_1_N"] = SelectionRange(5.0, 11.0);
	yAlignmentRangeX["R_1_N"] = SelectionRange(5.0, 11.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(4.0, 11.0);
}
