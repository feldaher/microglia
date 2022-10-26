
/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./OT_2D.h"

Cell_Definition immune_cell;
Cell_Definition neuron;
Cell_Definition skin_cell;
Cell_Definition neuropil_cell;



void create_immune_cell_type( void )
{
	immune_cell = cell_defaults;

	immune_cell.name = "immune_cell";
	immune_cell.type = immune_cell_ID;

	// turn off proliferation;

	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );

	immune_cell.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 0.0;

	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );

	// reduce o2 uptake

	immune_cell.phenotype.secretion.uptake_rates[0] *=
		parameters.doubles("immune_o2_relative_uptake"); // 0.1;

	// set apoptosis to survive 10 days (on average)

	immune_cell.phenotype.death.rates[apoptosis_index] =
		parameters.doubles("immune_apoptosis_rate"); // 1.0 / (10.0 * 24.0 * 60.0 );

	// turn on motility;
	immune_cell.phenotype.motility.is_motile = true;
	immune_cell.phenotype.motility.persistence_time =
		parameters.doubles("immune_motility_persistence_time"); // 10.0;
	immune_cell.phenotype.motility.migration_speed =
		parameters.doubles("immune_migration_speed"); // 1;
	immune_cell.phenotype.motility.migration_bias =
		parameters.doubles("immune_migration_bias"); // 0.5;

	immune_cell.phenotype.mechanics.cell_cell_adhesion_strength *=
		parameters.doubles("immune_relative_adhesion"); // 0.0;
	immune_cell.phenotype.mechanics.cell_cell_repulsion_strength *=
		parameters.doubles("immune_relative_repulsion"); // 5.0;

	// set functions

	immune_cell.functions.update_phenotype = NULL;
//	neuron.functions.custom_cell_rule = neuron_rule;
//	neuron.functions.update_migration_bias = neuron_motility;


	return;
}

void create_neuron_cell_type( void )
{
	neuron = cell_defaults;

	neuron.name = "neuron";
	neuron.type = neuron_ID;

	// turn off proliferation;

	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );

	neuron.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 0.0;

	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );

	// reduce o2 uptake

	neuron.phenotype.secretion.uptake_rates[0] *=
		parameters.doubles("neuron_o2_relative_uptake"); // 0.1;

	// set apoptosis to survive 10 days (on average)

	neuron.phenotype.death.rates[apoptosis_index] =
		parameters.doubles("neuron_apoptosis_rate"); // 1.0 / (10.0 * 24.0 * 60.0 );

	// turn on motility;
	neuron.phenotype.motility.is_motile = true;
	neuron.phenotype.motility.persistence_time =
		parameters.doubles("neuron_motility_persistence_time"); // 10.0;
	neuron.phenotype.motility.migration_speed =
		parameters.doubles("neuron_migration_speed"); // 1;
	neuron.phenotype.motility.migration_bias =
		parameters.doubles("neuron_migration_bias"); // 0.5;

	neuron.phenotype.mechanics.cell_cell_adhesion_strength *=
		parameters.doubles("neuron_relative_adhesion"); // 0.0;
	neuron.phenotype.mechanics.cell_cell_repulsion_strength *=
		parameters.doubles("neuron_relative_repulsion"); // 5.0;

	// set functions

	neuron.functions.update_phenotype = NULL;
//	neuron.functions.custom_cell_rule = neuron_rule;
//	neuron.functions.update_migration_bias = neuron_motility;


	return;
}


void create_skin_cell_type( void )
{
	skin_cell = cell_defaults;

	skin_cell.name = "skin_cell";
	skin_cell.type = skin_cell_ID;

	// turn off proliferation;

	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );

	skin_cell.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 0.0;

	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );



	// set apoptosis to survive 10 days (on average)

	skin_cell.phenotype.death.rates[apoptosis_index] =
		parameters.doubles("skin_cell_apoptosis_rate"); // 1.0 / (10.0 * 24.0 * 60.0 );

	// turn on motility;
	skin_cell.phenotype.motility.is_motile = false;

	skin_cell.phenotype.mechanics.cell_cell_adhesion_strength *=
		parameters.doubles("skin_cell_relative_adhesion"); // 0.0;
	skin_cell.phenotype.mechanics.cell_cell_repulsion_strength *=
		parameters.doubles("skin_cell_relative_repulsion"); // 5.0;

	// set functions

	skin_cell.functions.update_phenotype = NULL;
//	neuron.functions.custom_cell_rule = neuron_rule;
//	neuron.functions.update_migration_bias = neuron_motility;


	return;
}




void create_neuropil_cell_type( void )
{
	neuropil_cell = cell_defaults;

	neuropil_cell.name = "neuropil_cell";
	neuropil_cell.type = neuropil_cell_ID;

	// turn off proliferation;

	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );

	neuropil_cell.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 0.0;

	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );



	// set apoptosis to survive 10 days (on average)

	neuropil_cell.phenotype.death.rates[apoptosis_index] =
		parameters.doubles("skin_cell_apoptosis_rate"); // 1.0 / (10.0 * 24.0 * 60.0 );

	// turn on motility;
	neuropil_cell.phenotype.motility.is_motile = false;

	neuropil_cell.phenotype.mechanics.cell_cell_adhesion_strength *=
		parameters.doubles("skin_cell_relative_adhesion"); // 0.0;
	neuropil_cell.phenotype.mechanics.cell_cell_repulsion_strength *=
		parameters.doubles("skin_cell_relative_repulsion"); // 5.0;

	// set functions

	neuropil_cell.functions.update_phenotype = NULL;
//	neuron.functions.custom_cell_rule = neuron_rule;
//	neuron.functions.update_migration_bias = neuron_motility;


	return;
}

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the
	// same initial histogram of oncoprotein, even if threading means
	// that future division and other events are still not identical
	// for all runs
	SeedRandom( parameters.ints("random_seed") );

	// housekeeping

	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );

	// turn the default cycle model to live,
	// so it's easier to turn off proliferation

	cell_defaults.phenotype.cycle.sync_to_cycle_model( live );

	// Make sure we're ready for 2D

	cell_defaults.functions.set_orientation = up_orientation;
	// cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = false; // true;


	// set to no motility for cancer cells
	cell_defaults.phenotype.motility.is_motile = false;

	// set default motility parameters

	cell_defaults.phenotype.motility.is_motile = true;
	cell_defaults.phenotype.motility.persistence_time =
		parameters.doubles("default_motility_persistence_time"); // 10.0;
	cell_defaults.phenotype.motility.migration_speed =
		parameters.doubles("default_migration_speed"); // 1;
	cell_defaults.phenotype.motility.migration_bias =
		parameters.doubles("default_migration_bias"); // 0.5;

	// use default proliferation and death
	// turn off proliferation;

	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );

	cell_defaults.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 0.0;

	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );

	cell_defaults.phenotype.death.rates[apoptosis_index] =0;
	//Change the mechanics: cell_cell_repulsion_strength
	//cell_defaults.phenotype.mechanics.cell_cell_repulsion_strength=0;
	//cell_defaults.phenotype.mechanics.set_relative_equilibrium_distance(3.0);

	// set default uptake and secretion
	// oxygen
	cell_defaults.phenotype.secretion.secretion_rates[0] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[0] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[0] = 38;

	// immunostimulatory
	cell_defaults.phenotype.secretion.saturation_densities[1] = 1;

	// set the default cell type to o2-based proliferation with the effect of the
	// on oncoprotein, and secretion of the immunostimulatory factor


	// add the extra bit of "attachment" mechanics
	cell_defaults.functions.custom_cell_rule = extra_elastic_attachment_mechanics;

	//-------> set adhesive force to 0
	//cell_defaults.phenotype.mechanics.cell_cell_adhesion_strength=0;

	cell_defaults.name = "cancer cell";
	cell_defaults.type = 0;

	// add custom data

	Parameter<double> paramD;

	cell_defaults.custom_data.add_variable( "oncoprotein" , "dimensionless", 1.0 );
	paramD = parameters.doubles[ "elastic_coefficient" ];
	cell_defaults.custom_data.add_variable( "elastic coefficient" , paramD.units, paramD.value );
		// "1/min" , 0.01 );  /* param */
	cell_defaults.custom_data.add_variable( "kill rate" , "1/min" , 0 ); // how often it tries to kill
	cell_defaults.custom_data.add_variable( "attachment lifetime" , "min" , 0 ); // how long it can stay attached
	cell_defaults.custom_data.add_variable( "attachment rate" , "1/min" ,0 ); // how long it wants to wander before attaching



	// create the immune cell type
//	create_immune_cell_type();

	create_neuron_cell_type();
	create_immune_cell_type();
	create_skin_cell_type();
	create_neuropil_cell_type();


	return;
}


void setup_microenvironment( void )
{
	// set domain parameters

/* now this is in XML
	default_microenvironment_options.X_range = {-1000, 1000};
	default_microenvironment_options.Y_range = {-1000, 1000};
	default_microenvironment_options.simulate_2D = true;
*/

	// make sure to override and go back to 2D
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl;
		default_microenvironment_options.simulate_2D = true;
	}



	initialize_microenvironment();

	// Microenvironment thisMicroenvironment=NULL;
	//
	// thisMicroenvironment= get_default_microenvironment();


	// Set mechanics voxel size.
	double mechanics_voxel_size = 200;
	// Assume microenvironment is defined above somewhere.
	// Set up the PhysiCell mechanics data structure.
	Cell_Container* cell_container = create_cell_container_for_microenvironment(microenvironment, mechanics_voxel_size );


	return;
}

/************* TISSUE SETUP *************************/
void setup_tissue( void )
{
	// place a cluster of tumor cells at the center
//	cell_defaults.phenotype.geometry.radius=4;
//	double cell_radius = cell_defaults.phenotype.geometry.radius;
//	cell_defaults.phenotype.volume.total=12000;
	double cell_radius = std::cbrt(cell_defaults.phenotype.volume.total*0.2387);
//	double cell_spacing = 0.95 * 2.0 * cell_radius;

	std::cout << "cell radius" << cell_radius << std::endl;

	double PVZ_size =		parameters.doubles("PVZ_size"); // 250.0;
	double Patch_size =		parameters.doubles("Patch_size"); // 250.0;


	Cell* pCell = NULL;

//		std::vector<std::vector<double>> positions = create_cell_sphere_positions(cell_radius,PVZ_size);
//	std::vector<std::vector<double>> positions = create_cell_PVZ_positions(cell_radius,PVZ_size);
	std::vector<std::vector<double>> positions = create_cell_positions_from_file("PVZ.dat");
	std::cout << "creating " << positions.size() << " closely-packed neurons cells ... " << std::endl;

//	std::vector<std::vector<double>> positions_Immune = create_cell_NP_positions(cell_radius,Patch_size);
	std::vector<std::vector<double>> positions_Immune = create_cell_positions_from_file("microglia.dat");
	std::cout << "creating " << positions_Immune.size() << " immune cells patch... " << std::endl;




//	std::vector<std::vector<double>> positions_Skin = create_cell_skin_positions(cell_radius);
	std::vector<std::vector<double>> positions_Skin = create_cell_positions_from_file("skin.dat");
	std::cout << "creating " << positions_Skin.size() << " skin cells... " << std::endl;


	//	std::vector<std::vector<double>> positions_Skin = create_cell_skin_positions(cell_radius);
	std::vector<std::vector<double>> positions_RGC = create_cell_positions_from_file("neuropil.dat");
	std::cout << "creating " << positions_RGC.size() << " RGC cells... " << std::endl;

	

	int x_min=0;
	int x_max=0;
	int y_min=0;
	int y_max=0;


	x_min=default_microenvironment_options.X_range[0];
	x_max=PVZ_size;
	y_min=default_microenvironment_options.Y_range[0];
	y_max=default_microenvironment_options.Y_range[1];

	//SET NEURON POSITIONS
	for( int i=0; i < positions.size(); i++ )
	{
//		std::cout<<"creating neurons...."<<std::endl;
		pCell = create_cell(neuron); // tumor cell
		pCell->assign_position( positions[i] );
//		pCell->functions.add_cell_basement_membrane_interactions=true;

		//Set tissue fixed boundaries
		//at X min
		if(positions[i][0]<x_min+2*cell_radius){
			pCell->is_movable=false;
		}
		//at X max
		// if(positions[i][0]>x_max-2*cell_radius){
		// 	pCell->is_movable=false;
		// }
		// //at Y min
		if(positions[i][1]<y_min+2*cell_radius){
			pCell->is_movable=false;
		}

		//at Y max
		if(positions[i][1]>y_max-2*cell_radius){
			pCell->is_movable=false;
		}


	}



//SET SKIN CELLS POSITIONS
		for( int i=0; i < positions_Skin.size(); i++ )
		{
	//		std::cout<<"creating skin cells...."<<std::endl;
			pCell = create_cell(skin_cell); // tumor cell
			pCell->assign_position( positions_Skin[i] );
		}


//SET RGC CELLS POSITIONS
		for( int i=0; i < positions_RGC.size(); i++ )
		{
	//		std::cout<<"creating RGC...."<<std::endl;
			pCell = create_cell(neuropil_cell); // tumor cell
			pCell->assign_position( positions_RGC[i] );
		}

		
//SET IMMUNE CELLS POSITIONS
	for( int i=0; i < positions_Immune.size(); i++ )
	{
//		std::cout<<"creating immune cells...."<<std::endl;
		pCell = create_cell(immune_cell); // tumor cell
		pCell->assign_position( positions_Immune[i] );
	}



//SET ADHESION BETWEEN IMMUNE CELLS AND OTHER CELLS

	for( int i=0; i < all_cells->size() ; i++ )
	{
			if((*all_cells)[i]->type==immune_cell_ID){
				//for j!= i, attach
				for(int j=0;j < all_cells->size() ; j++ ){
					if(j!=i){
						attach_cells((*all_cells)[j],(*all_cells)[i]);
					}
				}
			}
	}

//SET ADHESION OF RGC BETWEEN THEM

	for( int i=0; i < all_cells->size() ; i++ )
	{
			if((*all_cells)[i]->type==neuropil_cell_ID){
				//for j!= i, attach
				for(int j=0;j < all_cells->size() ; j++ ){
					if((j!=i) & ((*all_cells)[j]->type==neuropil_cell_ID)){
						attach_cells((*all_cells)[j],(*all_cells)[i]);
					}
				}
			}
	}

//SET ADHESION BETWEEN RGC CELLS AND neurons

		for( int i=0; i < all_cells->size() ; i++ )
		{
				if((*all_cells)[i]->type==neuropil_cell_ID){
					//for j!= i, attach
					for(int j=0;j < all_cells->size() ; j++ ){
						if(j!=i){
							if((*all_cells)[j]->type==neuron_ID){
								attach_cells((*all_cells)[j],(*all_cells)[i]);
							}
						}
					}
				}
		}
	
	return;
	
}
/************* TISSUE GEOMETRY *************************/


/************* CELL POSITION ARRANGEMENT *************************/
std::vector<std::vector<double>> create_cell_PVZ_positions(double cell_radius, double PVZ_size)
{
	std::vector<std::vector<double>> positions;
	int xc=0,yc=0;
	double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*2;

	double x_min=0.0,x_max=0.0, y_min=0.0,y_max=0.0;
//	bool	create_inj=true;
	Parameter<bool> create_inj;
	create_inj = parameters.bools["create_injury"];

	std::vector<double> tempPoint(3,0.0);


	x_min=default_microenvironment_options.X_range[0];
	x_max=x_min+PVZ_size;
	y_min=default_microenvironment_options.Y_range[0];
	y_max=default_microenvironment_options.Y_range[1];

	// define the geometry of the hole
	double Injury_size =parameters.doubles("Injury_size"); // 250.0;
	double x_min_inj=parameters.doubles("x_min_inj");
	double y_min_inj=parameters.doubles("y_min_inj");

	double x_max_inj=x_min_inj+Injury_size;
	double y_max_inj=y_min_inj+Injury_size;

	std::cout << "x_min_inj: " << x_min_inj<< '\n';
	std::cout << "x_max_inj: " << x_max_inj<< '\n';
	std::cout << "y_min_inj: " << y_min_inj<< '\n';
	std::cout << "y_max_inj: " << y_max_inj<< '\n';


		for(double x=x_min;x<x_max;x+=x_spacing, xc++)
		{
			for(double y=y_min;y<y_max;y+=y_spacing, yc++)
			{
				tempPoint[0]=x + (yc%2) * 0.5 * cell_radius;
				tempPoint[1]=y;

				//Creating a tissue injury
				if (create_inj.value==true){
					if(tempPoint[0]<x_min_inj||tempPoint[0]>x_max_inj||tempPoint[1]<y_min_inj||tempPoint[1]>y_max_inj){
						positions.push_back(tempPoint);
							//std::cout << "Creating a cell in injured tissue" << '\n';
					}
				}else{
					positions.push_back(tempPoint);
					// std::cout << "Creating a cell in intact tissue" << '\n';
				}
			}
		}

	return positions;

}

std::vector<std::vector<double>> create_cell_positions_from_file(std::string filename)
{

	// READ FROM FILE THE POSITIONS AND CELL TYPE
	std::vector<std::vector<double>> positions;
	std::vector<double> tempPoint(2,0.0); // used to get x and y values from file

	std::ifstream in_data;
//	std::string filename;
//	filename="cells.dat";

	in_data.open (filename, std::ifstream::in);

	if(!in_data.is_open()) throw std::runtime_error("Could not open file");


	std::cout<<"file opened"<<std::endl;
	std::vector<std::string> row;
	std::string line, word, temp;
//	char* line;
	// char c = in_data.get();
	//
  // while (in_data.good()) {
  //   std::cout << c;
  //   c = in_data.get();
  // }

	if(in_data.good()){
		std::getline(in_data, line);// skip the column names in the first row
	//	std::cout << line;

		while (std::getline(in_data, line)) {
			row.clear();

	    // read an entire row and
	    // store it in a string variable 'line'
		//	std::cout << line;

			// used for breaking words
	    std::stringstream s(line);

	    // read every column data of a row and
	  	// store it in a string variable, 'word'
	    while (std::getline(s, word, ' ')) {

	        // add all the column data
	        // of a row to a vector
	        row.push_back(word);
					// tempPoint[]
			}
			tempPoint[0]=stod(row[0]); // x position
			tempPoint[1]=stod(row[1]); // y positions
		//	std::cout<< tempPoint[0]<<tempPoint[1]<<std::endl;

			positions.push_back(tempPoint);
		}
	}

	in_data.close();


	return positions;
}


std::vector<std::vector<double>> create_cell_NP_positions(double cell_radius, double Patch_size)
{
	std::vector<std::vector<double>> positions;
	int xc=0,yc=0;
	double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*3;

	double x_min=0.0,x_max=0.0, y_min=0.0,y_max=0.0;


	std::vector<double> tempPoint(3,0.0);


	//x_min=default_microenvironment_options.X_range[1]*0.3-Patch_size/2; //average position of microglia from exp.
x_min=-300;
	x_max=x_min+Patch_size;
	y_min=-Patch_size/2;
	y_max=Patch_size/2;


		for(double x=x_min;x<x_max;x+=x_spacing, xc++)
		{
			for(double y=y_min;y<y_max;y+=y_spacing, yc++)
			{
				tempPoint[0]=x + (yc%2) * 0.5 * cell_radius;
				tempPoint[1]=y;

				positions.push_back(tempPoint);
			}
		}

	return positions;

}

std::vector<std::vector<double>> create_cell_skin_positions(double cell_radius)
{
	std::vector<std::vector<double>> positions;
	int xc=0,yc=0;
	double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*2;

	double x_min=0.0,x_max=0.0, y_min=0.0,y_max=0.0;
//	bool	create_inj=true;
	Parameter<bool> create_inj;
	create_inj = parameters.bools["create_injury"];

	std::vector<double> tempPoint(3,0.0);


	x_min=default_microenvironment_options.X_range[1]-cell_radius*30;//thickness of the skin cells layer
	x_max=default_microenvironment_options.X_range[1];
	y_min=default_microenvironment_options.Y_range[0];
	y_max=default_microenvironment_options.Y_range[1];


		for(double x=x_min;x<x_max;x+=x_spacing, xc++)
		{
			for(double y=y_min;y<y_max;y+=y_spacing, yc++)
			{
				tempPoint[0]=x + (yc%2) * 0.5 * cell_radius;
				tempPoint[1]=y;
				positions.push_back(tempPoint);
			}
		}

	return positions;

}


/************* EXTRA CELL MECHANICS *************************/

void add_elastic_velocity( Cell* pActingOn, Cell* pAttachedTo , double elastic_constant )
{
	if(pActingOn->type!=skin_cell_ID){ // We don't want skin cells to move due to attractive forces
		std::vector<double> displacement = pAttachedTo->position - pActingOn->position;
		axpy( &(pActingOn->velocity) , elastic_constant , displacement );
	}


	return;
}

void extra_elastic_attachment_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	for( int i=0; i < pCell->state.neighbors.size() ; i++ )
	{
		add_elastic_velocity( pCell, pCell->state.neighbors[i], pCell->custom_data["elastic coefficient"] );
	}

	return;
}

void attach_cells( Cell* pCell_1, Cell* pCell_2 )
{
	#pragma omp critical
	{

	bool already_attached = false;
	for( int i=0 ; i < pCell_1->state.neighbors.size() ; i++ )
	{
		if( pCell_1->state.neighbors[i] == pCell_2 )
		{ already_attached = true; }
	}
	if( already_attached == false )
	{ pCell_1->state.neighbors.push_back( pCell_2 ); }

	already_attached = false;
	for( int i=0 ; i < pCell_2->state.neighbors.size() ; i++ )
	{
		if( pCell_2->state.neighbors[i] == pCell_1 )
		{ already_attached = true; }
	}
	if( already_attached == false )
	{ pCell_2->state.neighbors.push_back( pCell_1 ); }

	}

	return;
}

Cell* immune_cell_check_neighbors_for_attachment( Cell* pAttacker , double dt )
{
	std::vector<Cell*> nearby = pAttacker->cells_in_my_container();
	int i = 0;
	while( i < nearby.size() )
	{
		// don't try to kill yourself
		if( nearby[i] != pAttacker )
		{
			if( immune_cell_attempt_attachment( pAttacker, nearby[i] , dt ) )
			{ return nearby[i]; }
		}
		i++;
	}

	return NULL;
}

bool immune_cell_attempt_attachment( Cell* pAttacker, Cell* pTarget , double dt )
{

			attach_cells( pAttacker, pTarget );

	return true;
}



void immune_cell_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
	// if I'm docked
	if( pCell->state.neighbors.size() > 0 )
	{
		extra_elastic_attachment_mechanics( pCell, phenotype, dt );
	}

	immune_cell_check_neighbors_for_attachment( pCell , dt);

	return;
}




std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with flow cytometry coloring

	std::vector<std::string> output = false_cell_coloring_cytometry(pCell);
	std::string color = "black";

	if( pCell->phenotype.death.dead == true )
	{ return output; }

if(pCell->type == neuron_ID)
{ color = "red"; }
else if(pCell->type == immune_cell_ID)
{ color = "blue"; }
else if(pCell->type== skin_cell_ID)
{color = "black";}
else if(pCell->type== neuropil_cell_ID)
{color = "orange";}



output[0] = color;
output[2] = color;


	return output;
}
