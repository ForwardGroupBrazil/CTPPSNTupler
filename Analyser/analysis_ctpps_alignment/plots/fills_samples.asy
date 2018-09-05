struct Dataset
{
	string tag;
	int idx;
}

//----------------------------------------------------------------------------------------------------

struct FillData
{
	int fill;
	Dataset datasets[];
};

//----------------------------------------------------------------------------------------------------

FillData fill_data[];

void AddDataSet(string p, int idx = 0)
{
	int fill = (int) substr(p, find(p, "fill_")+5, 4);

	bool found = false;
	for (FillData fd : fill_data)
	{
		if (fd.fill == fill)
		{
			found = true;
			Dataset ds;
			ds.tag = p;
			ds.idx = idx;
			fd.datasets.push(ds);
		}
	}

	if (!found)
	{
		FillData fd;
		fd.fill = fill;
		Dataset ds;
		ds.tag = p;
		ds.idx = idx;
		fd.datasets.push(ds);

		fill_data.push(fd);
	}
}

//----------------------------------------------------------------------------------------------------

void InitDataSets(string selection = "period1")
{
	bool period1_ps = false;
	bool period1_ss = false;
	bool period2_ps = false;
	bool period2_ss = false;

	for (string bit : split(selection, ","))
	{
		if (bit == "period1_ps" || bit == "period1") period1_ps = true;
		if (bit == "period1_ss" || bit == "period1") period1_ss = true;
		if (bit == "period2_ps" || bit == "period2") period2_ps = true;
		if (bit == "period2_ss" || bit == "period2") period2_ss = true;
	}

	// period1, primary sample
	if (period1_ps)
	{
		AddDataSet("period1_physics_margin/fill_4947");
		AddDataSet("period1_physics_margin/fill_4953");
		AddDataSet("period1_physics_margin/fill_4961");
		AddDataSet("period1_physics_margin/fill_4964");
		
		AddDataSet("period1_physics/fill_4964");
		AddDataSet("period1_physics_margin/fill_4976");
		
		AddDataSet("period1_physics/fill_4985");
		AddDataSet("period1_physics/fill_4988");
		AddDataSet("period1_physics/fill_4990");
		AddDataSet("period1_physics/fill_5005");
		AddDataSet("period1_physics/fill_5013");
		AddDataSet("period1_physics/fill_5017");
		AddDataSet("period1_physics/fill_5020");
		AddDataSet("period1_physics/fill_5021");
		AddDataSet("period1_physics/fill_5024");
		AddDataSet("period1_physics/fill_5026");
		AddDataSet("period1_physics/fill_5027");
		AddDataSet("period1_physics/fill_5028");
		AddDataSet("period1_physics/fill_5029");
		AddDataSet("period1_physics/fill_5030");
		AddDataSet("period1_physics/fill_5038");
		AddDataSet("period1_physics/fill_5043");
		AddDataSet("period1_physics/fill_5045");
		AddDataSet("period1_physics/fill_5048");
		AddDataSet("period1_physics/fill_5052");
		
		AddDataSet("period1_physics/fill_5261");
		AddDataSet("period1_physics/fill_5264");
		AddDataSet("period1_physics/fill_5265");
		AddDataSet("period1_physics/fill_5266");
		AddDataSet("period1_physics/fill_5267");
		AddDataSet("period1_physics/fill_5274");
		AddDataSet("period1_physics/fill_5275");
		AddDataSet("period1_physics/fill_5276");
		AddDataSet("period1_physics/fill_5277");
		AddDataSet("period1_physics/fill_5279");
		AddDataSet("period1_physics/fill_5287");
		AddDataSet("period1_physics/fill_5288");
	}
	
	// period1, secondary sample
	if (period1_ss)
	{
		AddDataSet("period1_physics_margin_ss/fill_4947", 1);
		AddDataSet("period1_physics_margin_ss/fill_4961", 1);
		AddDataSet("period1_physics_margin_ss/fill_4976", 1);
		
		AddDataSet("period1_physics_ss/fill_4985", 1);
		AddDataSet("period1_physics_ss/fill_4988", 1);
		AddDataSet("period1_physics_ss/fill_4990", 1);
		AddDataSet("period1_physics_ss/fill_5005", 1);
		AddDataSet("period1_physics_ss/fill_5013", 1);
		AddDataSet("period1_physics_ss/fill_5017", 1);
		AddDataSet("period1_physics_ss/fill_5020", 1);
		AddDataSet("period1_physics_ss/fill_5021", 1);
		AddDataSet("period1_physics_ss/fill_5024", 1);
		AddDataSet("period1_physics_ss/fill_5026", 1);
		AddDataSet("period1_physics_ss/fill_5027", 1);
		AddDataSet("period1_physics_ss/fill_5028", 1);
		AddDataSet("period1_physics_ss/fill_5029", 1);
		AddDataSet("period1_physics_ss/fill_5030", 1);
		AddDataSet("period1_physics_ss/fill_5038", 1);
		AddDataSet("period1_physics_ss/fill_5043", 1);
		AddDataSet("period1_physics_ss/fill_5045", 1);
		AddDataSet("period1_physics_ss/fill_5048", 1);
		AddDataSet("period1_physics_ss/fill_5052", 1);
		AddDataSet("period1_physics_ss/fill_5261", 1);
		AddDataSet("period1_physics_ss/fill_5267", 1);
		AddDataSet("period1_physics_ss/fill_5277", 1);
		AddDataSet("period1_physics_ss/fill_5279", 1);
		AddDataSet("period1_physics_ss/fill_5287", 1);
		AddDataSet("period1_physics_ss/fill_5288", 1);
	}

	// period2, primary sample
	if (period2_ps)
	{
		AddDataSet("period2_physics/fill_5393");
		AddDataSet("period2_physics/fill_5401");
		AddDataSet("period2_physics/fill_5405");
		AddDataSet("period2_physics/fill_5406");
		AddDataSet("period2_physics/fill_5418");
		AddDataSet("period2_physics/fill_5421");
		AddDataSet("period2_physics/fill_5423");
		AddDataSet("period2_physics/fill_5424");
		AddDataSet("period2_physics/fill_5427");
		AddDataSet("period2_physics/fill_5433");
		AddDataSet("period2_physics/fill_5437");
		AddDataSet("period2_physics/fill_5439");
		AddDataSet("period2_physics/fill_5441");
		AddDataSet("period2_physics/fill_5442");
		AddDataSet("period2_physics/fill_5443");
		AddDataSet("period2_physics/fill_5446");
		AddDataSet("period2_physics/fill_5448");
		AddDataSet("period2_physics/fill_5450");
		AddDataSet("period2_physics/fill_5451");
	}
}
