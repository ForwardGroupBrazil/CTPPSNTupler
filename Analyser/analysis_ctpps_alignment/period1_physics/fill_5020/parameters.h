void ApplySettings()
{
	lhc_fill = 5020;

	cut1_a = -1.008; cut1_c = -1.343; cut1_si = 0.20;
	cut2_a = -0.978; cut2_c = +0.544; cut2_si = 0.20;
	cut3_a = -1.166; cut3_c = +0.772; cut3_si = 0.15;
	cut4_a = -1.145; cut4_c = +0.037; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(11.0, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(9.3, 14.5);
	selectionRangesX["R_1_N"] = SelectionRange(8.9, 15.5);
	selectionRangesX["R_1_F"] = SelectionRange(8.5, 14.5);

	yAlignmentRangeX["L_1_F"] = SelectionRange(6.5, 12.0);
	yAlignmentRangeX["L_1_N"] = SelectionRange(6.5, 12.0);
	yAlignmentRangeX["R_1_N"] = SelectionRange(5.0, 10.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(4.0, 11.0);
}
