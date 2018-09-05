void ApplySettings()
{
	lhc_fill = 5038;

	cut1_a = -1.006; cut1_c = -1.009; cut1_si = 0.20;
	cut2_a = -0.978; cut2_c = +0.532; cut2_si = 0.20;
	cut3_a = -1.163; cut3_c = +0.507; cut3_si = 0.15;
	cut4_a = -1.141; cut4_c = +0.042; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(11.3, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(10.2, 14.5);
	selectionRangesX["R_1_N"] = SelectionRange(8.8, 15.5);
	selectionRangesX["R_1_F"] = SelectionRange(8.2, 14.5);

	yAlignmentRangeX["L_1_F"] = SelectionRange(7.5, 11.0);
	yAlignmentRangeX["L_1_N"] = SelectionRange(7.0, 11.0);
	yAlignmentRangeX["R_1_N"] = SelectionRange(5.0, 11.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(4.0, 11.0);
}
