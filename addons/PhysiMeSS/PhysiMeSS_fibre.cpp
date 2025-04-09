#include "PhysiMeSS_fibre.h"
#include "PhysiMeSS_cell.h"
#include <algorithm>

bool isFibre(PhysiCell::Cell* pCell) 
{
    const auto agentname = std::string(pCell->type_name);
    const auto ecm = std::string("ecm");
    const auto matrix = std::string("matrix");
    const auto fiber = std::string("fiber");
    const auto fibre = std::string("fibre");
    const auto rod = std::string("rod");

    return (agentname.find(ecm) != std::string::npos ||
        agentname.find(matrix) != std::string::npos ||
        agentname.find(fiber) != std::string::npos ||
        agentname.find(fibre) != std::string::npos ||
        agentname.find(rod) != std::string::npos
    );
}

bool isFibre(PhysiCell::Cell_Definition * cellDef)
{
    const auto agentname = std::string(cellDef->name);
    const auto ecm = std::string("ecm");
    const auto matrix = std::string("matrix");
    const auto fiber = std::string("fiber");
    const auto fibre = std::string("fibre");
    const auto rod = std::string("rod");

    return (agentname.find(ecm) != std::string::npos ||
        agentname.find(matrix) != std::string::npos ||
        agentname.find(fiber) != std::string::npos ||
        agentname.find(fibre) != std::string::npos ||
        agentname.find(rod) != std::string::npos
    );
}

std::vector<PhysiCell::Cell_Definition*>* getFibreCellDefinitions() {
    std::vector<PhysiCell::Cell_Definition*>* result = new std::vector<PhysiCell::Cell_Definition*>();
    PhysiCell::Cell_Definition* pCD;
    
    
    for (auto& cd_name: PhysiCell::cell_definitions_by_name) {
        if (isFibre(cd_name.second)) {
            result->push_back(cd_name.second);
        }
    }
    
    return result;
}

PhysiMeSS_Fibre::PhysiMeSS_Fibre() 
{
    fibres_crosslinkers.clear();
    fibres_crosslink_point.clear();
    X_crosslink_count = 0;
    fail_count = 0;
}

void PhysiMeSS_Fibre::assign_fibre_orientation() 
{ 
    mLength = PhysiCell::NormalRandom(this->custom_data["fibre_length"], this->custom_data["length_normdist_sd"]) / 2.0;//mLength initialization
    mRadius = this->custom_data["fibre_radius"];
    this->assign_orientation();
    if (default_microenvironment_options.simulate_2D) {
        if (this->custom_data["anisotropic_fibres"] > 0.5){
            double theta = PhysiCell::NormalRandom(this->custom_data["fibre_angle_theta"], this->custom_data["angle_normdist_sd"]);
            this->state.orientation[0] = cos(theta);
            this->state.orientation[1] = sin(theta);
        }
        else{
            this->state.orientation = PhysiCell::UniformOnUnitCircle();
        }
        this->state.orientation[2] = 0.0;
    }
    else {
        if (this->custom_data["anisotropic_fibres"] > 0.5) {
            double theta = PhysiCell::NormalRandom(this->custom_data["fibre_angle_theta"], this->custom_data["angle_normdist_sd"]);
            double phi = PhysiCell::NormalRandom(this->custom_data["fibre_angle_phi"], this->custom_data["angle_normdist_sd"]);
            this->state.orientation[0] = cos(theta)*sin(phi);
            this->state.orientation[1] = sin(theta)*sin(phi);
            this->state.orientation[2] = cos(phi);
        }
        else {
            this->state.orientation = PhysiCell::UniformOnUnitSphere();
        }
        
    }
    //###########################################//
    //   this bit a hack for PacMan and maze	 //
    //###########################################//
    if (this->type_name == "fibre_vertical") {
        this->state.orientation[0] = 0.0;
        this->state.orientation[1] = 1.0;
        this->state.orientation[2] = 0.0;
    }
    if (this->type_name == "fibre_horizontal") {
        this->state.orientation[0] = 1.0;
        this->state.orientation[1] = 0.0;
        this->state.orientation[2] = 0.0;
    }
    //###########################################// 
}

void PhysiMeSS_Fibre::check_out_of_bounds(std::vector<double>& position)
{
    double Xmin = BioFVM::get_default_microenvironment()->mesh.bounding_box[0]; 
	double Ymin = BioFVM::get_default_microenvironment()->mesh.bounding_box[1]; 
	double Zmin = BioFVM::get_default_microenvironment()->mesh.bounding_box[2]; 

	double Xmax = BioFVM::get_default_microenvironment()->mesh.bounding_box[3]; 
	double Ymax = BioFVM::get_default_microenvironment()->mesh.bounding_box[4]; 
	double Zmax = BioFVM::get_default_microenvironment()->mesh.bounding_box[5]; 
    
    // start and end points of a fibre are calculated from fibre center
    double xs = position[0] - this->mLength * this->state.orientation[0];
    double xe = position[0] + this->mLength * this->state.orientation[0];
    double ys = position[1] - this->mLength * this->state.orientation[1];
    double ye = position[1] + this->mLength * this->state.orientation[1];
    double zs = position[2] - this->mLength * this->state.orientation[2];
    double ze = position[2] + this->mLength * this->state.orientation[2];
   
    /* check whether a fibre end point leaves the domain and if so initialise fibre again
                assume user placed the centre of fibre within the domain so reinitialise orientation,
                break after 10 failures
     */
    
    if (default_microenvironment_options.simulate_2D) {
        zs = 0;
        ze = 0;
    }
    // Anisotropic fibres have similar orientation, therefore if one of the endpoint is out of bound, no fibres will be set at this position
    if (this->custom_data["anisotropic_fibres"]) {
        if (xs < Xmin || xe > Xmax || xe < Xmin || xs > Xmax ||
            ys < Ymin || ye > Ymax || ye < Ymin || ys > Ymax ||
            zs < Zmin || ze > Zmax || ze < Zmin || zs > Zmax ) {
            fail_count = 10;
        }
        
    }
    //Otherwise, if fibre orientation is not constrained an other random orientation is set if out of bound
    else{
        while (fail_count < 10) {
            if (xs < Xmin || xe > Xmax || xe < Xmin || xs > Xmax ||
                ys < Ymin || ye > Ymax || ye < Ymin || ys > Ymax ||
                zs < Zmin || ze > Zmax || ze < Zmin || zs > Zmax ) {
                    fail_count++;

                    this->state.orientation = PhysiCell::UniformOnUnitSphere();
                    
                    xs = position[0] - mLength * this->state.orientation[0];
                    xe = position[0] + mLength * this->state.orientation[0];
                    ys = position[1] - mLength * this->state.orientation[1];
                    ye = position[1] + mLength * this->state.orientation[1];
                    zs = position[2] - mLength * this->state.orientation[2];
                    ze = position[2] + mLength * this->state.orientation[2];

                    if (default_microenvironment_options.simulate_2D) {
                        this->state.orientation[2] = 0;
                        zs = 0;
                        ze = 0;
                        
                    }
            }
            else {
                break;
            }
        }
    }
}


//Modification of terminology to prevent confusion in the type of the object
void PhysiMeSS_Fibre::add_potentials_from_cell(PhysiMeSS_Cell* cell) {

    if ( !(this->phenotype.motility.is_motile) || this->X_crosslink_count >= 2) {
        return;
    }

    nearest_point_on_fibre(cell->position, this->displacement);
    if (default_microenvironment_options.simulate_2D) {
        add_potentials_from_cell_2D(cell);
    }

    else if (!default_microenvironment_options.simulate_2D) {
        add_potentials_from_cell_3D(cell);
    }
}


void PhysiMeSS_Fibre::add_potentials_from_cell_2D(PhysiMeSS_Cell* cell) {
    double distance = 0.0;
    this->displacement[2] = 0.0;

    distance = std::max(sqrt(dot_product(this->displacement,this->displacement)), 0.00001);

    double R = cell->phenotype.geometry.radius + mRadius;
    if (distance > R) {
        return;
    }     
    else {
        //computation point collision fibre cell
        std::vector<double> point_of_impact(3, 0.0); 
        std::transform((*cell).position.begin(), (*cell).position.begin() + 2,
                        this->displacement.begin(),
                        point_of_impact.begin(),
                        [](double pos, double disp) { return pos - disp; });

        //center of rotation
        std::vector<double> center(3,0.0); 
        if ( X_crosslink_count == 0 ) {
            std::copy(this->position.begin(),
                    this->position.end(),
                    center.begin()); 
        }
        else if ( X_crosslink_count == 1) {
            std::copy(this->fibres_crosslink_point.begin(),
                this->fibres_crosslink_point.end(),
                center.begin()); 
        }

        std::vector<double> center_to_impact_vector(3,0.0);
        std::transform(center.begin(), center.begin() + 2,
                        point_of_impact.begin(),
                        center_to_impact_vector.begin(),
                        [](double pos, double imp) { return imp - pos; });
        double norm_center_to_impact = sqrt(dot_product(center_to_impact_vector, center_to_impact_vector));
        
        std::vector<double> old_orientation(3,0.0);
        std::transform(center_to_impact_vector.begin(), center_to_impact_vector.begin() + 2,
                        old_orientation.begin(),
                        [norm_center_to_impact](double pos) { return pos/norm_center_to_impact; });

        if (cell->custom_data["fibre_rotation"] > 0.5) {
            if (norm_center_to_impact < 1e-10) { 
                    return; 
            }

            std::vector<double> torque(3,0.0);
            torque = cross_product(center_to_impact_vector,(*cell).phenotype.motility.motility_vector);

            // moment_arm_magnitude = length of pivot if orientation is a unit vector
            double torque_magnitude = sqrt(dot_product(torque,torque));
            double fibre_length = 2 * mLength;     
            double moment_of_inertia = fibre_length * fibre_length / 12; //Formula for a rod
                    
            double angular_acceleration = custom_data["fibre_sticky"] * torque_magnitude / moment_of_inertia;
            double angle = angular_acceleration;
            
            if( angular_acceleration < 1e-16 ) {
                    return;
            }    
            //torque normalisation 
            std::transform(torque.begin(), torque.end(),
                        torque.begin(),
                        [torque_magnitude](double tor) { return tor/torque_magnitude; });
                    
            double c = std::cos(angle);
            double s = std::sin(angle);

            std::vector<double> v_cross_k(3,0.0);
            // (k x v)
            v_cross_k = cross_product(torque,old_orientation);
            std::vector<double> new_orientation(3,0.0);
                    
            std::transform(old_orientation.begin(), old_orientation.begin() + 2,
                            v_cross_k.begin(),
                            new_orientation.begin(),
                            [c, s](double or_val, double v_val) { return or_val * c + v_val * s; });

            state.orientation = new_orientation;

            normalize(&state.orientation);

            //CASE WHERE THE CROSSLINK POINT DIFFERENT FROM CENTER IS THE CENTER OF ROTATION
            double distance_center_position = PhysiCell::dist(this->position, center);
            
            std::vector<double> center_to_position_vector(3,0.0);
            std::transform(this->position.begin(), this->position.begin() + 2,
                            center.begin(),
                            center_to_position_vector.begin(),
                            [](double p, double c) { return p - c; });

            if ( distance_center_position > 1e-10) {
                if (dot_product(state.orientation,center_to_position_vector)/distance_center_position > 0) {
                    for (int i=0; i<2; i++ ) {
                        this-> position[i] = state.orientation[i] * distance_center_position + center[i];
                    }
                }
                else { 
                    for (int i=0; i<2; i++ ) {
                        this-> position[i] = -state.orientation[i]* distance_center_position + center[i];
                    }
                }
            }   
        }
        //cell-fibre pushing only if fibre has no crosslinks or if it has crosslink and its crosslinkers has no more than 2
        if ( X_crosslink_count == 0 || 
            ((this->fibres_crosslinkers).size() == 1 && static_cast<PhysiMeSS_Fibre*>((this->fibres_crosslinkers[0]))->fibres_crosslinkers.size() < 3)) {       
            if (cell->custom_data["fibre_pushing"] > 0.5) {
                // as per PhysiCell
                static double simple_pressure_scale = 0.027288820670331;//what does this do precisely?
                double temp_r = 1.0 - distance/R;
                temp_r *= temp_r;
                state.simple_pressure += (temp_r / simple_pressure_scale);
                
                double effective_repulsion = sqrt(phenotype.mechanics.cell_cell_repulsion_strength *
                                            (*cell).phenotype.mechanics.cell_cell_repulsion_strength);
                temp_r *= effective_repulsion;
                
                if (fabs(temp_r) > 1e-16) {// if temp_r very small <=> almost no overlap between fibre and cell, the overlap must be "significant" for rotation to occur
                    temp_r /= distance; //temp_r = (1.0 - distance/R)²/distance      
                    naxpy(&velocity, temp_r, displacement); // Performs the operation: this->velocity = this->velocity - temp_r * displacement
                }
            }
        }
    }
}

/*
In 3D simulation, this functions modifies the velocity and orientation of a fibre in function of the cells in contact with it.
*/
void PhysiMeSS_Fibre::add_potentials_from_cell_3D(PhysiMeSS_Cell* cell) {
    double distance = 0.0;
    distance = std::max(sqrt(dot_product(this->displacement,this->displacement)), 0.00001);
    double R = cell->phenotype.geometry.radius + this->mRadius;

    if (distance > R){
        return; // no interaction
    }
    else if (distance <= R) { 
        std::vector<double> point_of_impact(3, 0.0); // displacement vector is the vector from cell center to the nearest position on the fibre,
        // point of impact is the point on the fibre vector nearest to the cell
        std::transform((*cell).position.begin(), (*cell).position.end(),
                        this->displacement.begin(),
                        point_of_impact.begin(),
                        [](double pos, double disp) { return pos - disp; });
        
         //center of rotation
         std::vector<double> center(3,0.0); 
         if ( X_crosslink_count == 0 ) {
             std::copy(this->position.begin(),
                     this->position.end(),
                     center.begin()); 
         }
         else if ( X_crosslink_count == 1) {
             std::copy(this->fibres_crosslink_point.begin(),
                 this->fibres_crosslink_point.end(),
                 center.begin()); 
         }
         std::vector<double> center_to_impact_vector(3,0.0);
         std::transform(center.begin(), center.end(),
                         point_of_impact.begin(),
                         center_to_impact_vector.begin(),
                         [](double pos, double imp) { return imp - pos; });
         double norm_center_to_impact = sqrt(dot_product(center_to_impact_vector, center_to_impact_vector));
         
         std::vector<double> old_orientation(3,0.0);
         std::transform(center_to_impact_vector.begin(), center_to_impact_vector.end(),
                         old_orientation.begin(),
                         [norm_center_to_impact](double pos) { return pos/norm_center_to_impact; });
        
        if (cell->custom_data["fibre_rotation"] > 0.5) {
            if (norm_center_to_impact < 1e-10) { 
                return; 
            }

            //torque = orientation x motility
            std::vector<double> torque(3,0.0);
            torque = cross_product(center_to_impact_vector,(*cell).phenotype.motility.motility_vector);

            // moment_arm_magnitude = length of pivot if orientation is a unit vector
            double torque_magnitude = sqrt(dot_product(torque,torque));
            double fibre_length = 2 * mLength;

            double moment_of_inertia = fibre_length * fibre_length / 12; //Formula for a rod
            double angular_acceleration = custom_data["fibre_sticky"] * torque_magnitude / moment_of_inertia;
            double angle = angular_acceleration;
                
            if( angular_acceleration < 1e-16 ) {
                // if angular acceleration is small no tilt
                return;
            }

            //torque normalization
            std::transform(torque.begin(), torque.end(),
                        torque.begin(),
                        [torque_magnitude](double tor) { return tor/torque_magnitude; });
            

            double c = std::cos(angle);
            double s = std::sin(angle);
            
            // (k x v)
            std::vector<double> v_cross_k(3,0.0);
            v_cross_k = cross_product(torque,old_orientation);

            // Because pivot dot v=0 if pivot is perpendicular
            // the full Rodrigues formula simplifies (https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula):
            // v_new = vc + (k x v)s

            std::vector<double> new_orientation(3,0.0);
            std::transform(old_orientation.begin(), old_orientation.end(),
                            v_cross_k.begin(),
                            new_orientation.begin(),
                            [c, s](double or_val, double v_val) { return or_val * c + v_val * s; });

            state.orientation = new_orientation;
            normalize(&state.orientation);

            //CASE WHERE THE CROSSLINK POINT DIFFERENT FROM CENTER IS THE CENTER OF ROTATION
            double distance_center_position = PhysiCell::dist(this->position, center);
            
            std::vector<double> center_to_position_vector(3,0.0);
            std::transform(this->position.begin(), this->position.end(),
                            center.begin(),
                            center_to_position_vector.begin(),
                            [](double p, double c) { return p - c; });

            if ( distance_center_position > 1e-10) {
                if (dot_product(state.orientation,center_to_position_vector)/distance_center_position > 0) {
                    for (int i=0; i<3; i++ ) {
                        this-> position[i] = state.orientation[i] * distance_center_position + center[i];
                    }
                }
                else { 
                    for (int i=0; i<3; i++ ) {
                        this-> position[i] = -state.orientation[i]* distance_center_position + center[i];
                    }
                }
            }   
        }

        if ( X_crosslink_count == 0 || 
            ((this->fibres_crosslinkers).size() == 1 && static_cast<PhysiMeSS_Fibre*>((this->fibres_crosslinkers[0]))->fibres_crosslinkers.size() < 3)) {       
            if (cell->custom_data["fibre_pushing"] > 0.5) {
                // as per PhysiCell
                static double simple_pressure_scale = 0.027288820670331;//what does this do precisely?
                double temp_r = 1 - distance/R;
                temp_r *= temp_r;
                // add the relative pressure contribution NOT SURE IF NEEDED
                state.simple_pressure += (temp_r / simple_pressure_scale);

                double effective_repulsion = sqrt(phenotype.mechanics.cell_cell_repulsion_strength *
                (*cell).phenotype.mechanics.cell_cell_repulsion_strength);
                temp_r *= effective_repulsion;

                if (fabs(temp_r) > 1e-16) {// if temp_r very small <=> almost no overlap between fibre and cell, the overlap must be "significant" for rotation to occur
                    temp_r /= distance; //temp_r = (1.0 - distance/R)²/distance      
                    naxpy(&velocity, temp_r, displacement); // Performs the operation: this->velocity = this->velocity - temp_r * displacement
                }
            }
        }
    }
}


void PhysiMeSS_Fibre::add_potentials_from_fibre(PhysiMeSS_Fibre* other_fibre) 
{
    /* probably want something here to model tension along fibres
        * this will be strong tension along the fibre for those fibres with a crosslink
        * and weak tension with background ECM */
    return;
}

void PhysiMeSS_Fibre::register_fibre_voxels() {

    int voxel;
    int voxel_size = this->get_container()->underlying_mesh.dx; // note this must be the same as the mechanics_voxel_size
    int test = 2.0 * this->mLength / voxel_size; //allows us to sample along the fibre

    std::vector<double> fibre_start(3, 0.0);
    std::vector<double> fibre_end(3, 0.0);
    for (unsigned int i = 0; i < 3; i++) {
        fibre_start[i] = this->position[i] - this->mLength * this->state.orientation[i];
        fibre_end[i] = this->position[i] + this->mLength * this->state.orientation[i];
    }
    // first add the voxel of the fibre end point
    voxel = this->get_container()->underlying_mesh.nearest_voxel_index(fibre_end);
    //std::cout << voxel << " " ;
    physimess_voxels.push_back(voxel);
    if (std::find(this->get_container()->agent_grid[voxel].begin(),
                    this->get_container()->agent_grid[voxel].end(),
                    this) == this->get_container()->agent_grid[voxel].end()) {
        this->get_container()->agent_grid[voxel].push_back(this);
    }
    // then walk along the fibre from fibre start point sampling and adding voxels as we go
    std::vector<double> point_on_fibre(3, 0.0);
    for (unsigned int j = 0; j < test + 1; j++) {
        for (unsigned int i = 0; i < 3; i++) {
            point_on_fibre[i] = fibre_start[i] + j * voxel_size * this->state.orientation[i];
        }
        voxel = this->get_container()->underlying_mesh.nearest_voxel_index(point_on_fibre);
        //std::cout << voxel << " " ;
        physimess_voxels.push_back(voxel);
        if (std::find(this->get_container()->agent_grid[voxel].begin(),
                        this->get_container()->agent_grid[voxel].end(),
                        this) == this->get_container()->agent_grid[voxel].end()) {
            this->get_container()->agent_grid[voxel].push_back(this);
        }
    }
    //std::cout << std::endl;

    physimess_voxels.sort();
    physimess_voxels.unique();
}

void PhysiMeSS_Fibre::deregister_fibre_voxels() 
{
    int centre_voxel = this->get_container()->underlying_mesh.nearest_voxel_index(this->position);
    for (int voxel: physimess_voxels) {
        if (voxel != centre_voxel) {
            this->get_container()->remove_agent_from_voxel(this, voxel);
        }
    }
}


/*
This function modifies computes the point on the central axis of the fibre closest to the argument point
*/
std::vector<double> PhysiMeSS_Fibre::nearest_point_on_fibre(std::vector<double> point, std::vector<double> &displacement) 
{
    // don't bother if the "fibre_agent" is not a fibre
    if (!isFibre(this)) { return displacement; }

    double fibre_length = 2 * this->mLength;
    // vector pointing from one endpoint of "fibre_agent" to "point"
    std::vector<double> fibre_to_agent(3, 0.0);
    // scalar product fibre_to_agent * fibre_vector
    double fibre_to_agent_dot_endpoint_to_center = 0;
    std::vector<double> endpoint(3,0.0);
    double distance = 0;

    for (unsigned int i = 0; i < 3; i++) {
        endpoint[i] = this->position[i] - this->mLength * this->state.orientation[i];
        fibre_to_agent[i] = point[i] - endpoint[i];
        fibre_to_agent_dot_endpoint_to_center += fibre_to_agent[i] * 2 * (position[i] - endpoint[i]);//Magnitude of the vector fiber agent projected on the vector orientation of the fiber
    }
    // First 2 cases: one of the endpoint is the closest point on the fibre from the point
    // "point" is closest to the selected endpoint of this fibre
    if (fibre_to_agent_dot_endpoint_to_center < 0.) {
        for (int i = 0; i < 3; i++) {
            displacement[i] = fibre_to_agent[i];
        }

    }
        // “point” is closest to the other endpoint of this fibre
    else if (fibre_to_agent_dot_endpoint_to_center > fibre_length * fibre_length) {
        for (unsigned int i = 0; i < 3; i++) {
            displacement[i] = point[i] - (this->position[i]
                                            + this->mLength * this->state.orientation[i]);
        }
    }
       
       // Final case: “point” is closest to a point along this fibre
    else {
        for (unsigned int i = 0; i < 3; i++) {
            displacement[i] = fibre_to_agent[i] - (position[i] - endpoint[i])/mLength * fibre_to_agent_dot_endpoint_to_center/fibre_length;
        }
    }
    // the function returns the displacement vector
    return displacement;
}



void PhysiMeSS_Fibre::check_fibre_crosslinks(PhysiMeSS_Fibre *fibre_neighbor) {

    if (this == fibre_neighbor) { return; }

    if (isFibre(this) && isFibre(fibre_neighbor)) {

        // fibre endpoints
        std::vector<double> point1(3, 0.0);
        std::vector<double> point2(3, 0.0);
        std::vector<double> point3(3, 0.0);
        std::vector<double> point4(3, 0.0);
        for (int i = 0; i < 3; i++) {
            // endpoints of "this" fibre
            point1[i] = this->position[i] - mLength * this->state.orientation[i];
            point2[i] = this->position[i] + mLength * this->state.orientation[i];
            // endpoints of "neighbor" fibre
            point3[i] = fibre_neighbor->position[i] - fibre_neighbor->mLength * fibre_neighbor->state.orientation[i];
            point4[i] = fibre_neighbor->position[i] + fibre_neighbor->mLength * fibre_neighbor->state.orientation[i];
        }

        //vectors between fibre endpoints
        std::vector<double> p1_to_p2(3, 0.0);
        std::vector<double> p3_to_p4(3, 0.0);
        std::vector<double> p1_to_p3(3, 0.0);
        std::vector<double> centre_to_centre(3, 0.0);
        for (int i = 0; i < 3; i++) {
            // "fibre" fibre vector
            p1_to_p2[i] = point2[i] - point1[i];
            // "neighbor" fibre vector
            p3_to_p4[i] = point4[i] - point3[i];
            // vector from "fibre" to "neighbor"
            p1_to_p3[i] = point3[i] - point1[i];
            // vector between fibre centres
            centre_to_centre[i] = fibre_neighbor->position[i] - this->position[i];
        }

        double co_radius = this->mRadius + fibre_neighbor->mRadius;
        double co_length = this->mLength + fibre_neighbor->mLength;
        std::vector<double> zero(3, 0.0);
        double distance = PhysiCell::dist(zero, centre_to_centre);
        normalize(&centre_to_centre);

        /* test if fibres intersect
           (1) if fibres are coplanar and parallel:
           the cross product of the two fibre vectors is zero
           [(P2 - P1) x (P4 - P3)].[(P2 - P1) x (P4 - P3)] = 0 */
        std::vector<double> FCP = cross_product(p1_to_p2, p3_to_p4);
        /*  coplanar parallel fibres could intersect if colinear
            i.e. the orientation of the fibres are parallel or
            antiparallel to the centre_to_centre vector and
            distance between fibre centres is less than their co_length */
        if (dot_product(FCP,FCP) == 0 &&
            (centre_to_centre == this->state.orientation ||
            centre_to_centre == -1.0 * this->state.orientation) &&
            distance <= co_length) {
            //std::cout << "fibre " << fibre->ID << " crosslinks with parallel colinear fibre " <<  (*fibre_neighbor).ID << std::endl;
            if (std::find(this->fibres_crosslinkers.begin(), this->fibres_crosslinkers.end(), (fibre_neighbor)) ==
                this->fibres_crosslinkers.end()) {
                this->fibres_crosslinkers.push_back(fibre_neighbor);
            }
            this->fibres_crosslink_point = this->position + this->mLength * centre_to_centre;
        }
        /* (2) parallel fibres may sit on top of one another
            we check the distance between fibre end points and
            the nearest point on neighbor fibre to see if they do */
        std::vector<double> displacement(3, 0.0);
        fibre_neighbor->nearest_point_on_fibre(point1, displacement);
        double test_point1 = PhysiCell::dist(zero, displacement);
        fibre_neighbor->nearest_point_on_fibre(point2, displacement);
        double test_point2 = PhysiCell::dist(zero, displacement);
        this->nearest_point_on_fibre(point3, displacement);
        double test_point3 = PhysiCell::dist(zero, displacement);
        this->nearest_point_on_fibre(point4, displacement);
        double test_point4 = PhysiCell::dist(zero, displacement);
        if (std::abs(test_point1) <= co_radius ||
            std::abs(test_point2) <= co_radius ||
            std::abs(test_point3) <= co_radius ||
            std::abs(test_point4) <= co_radius &&
            centre_to_centre != this->state.orientation &&
            centre_to_centre != -1.0 * this->state.orientation) {
            //std::cout << "fibre " << fibre->ID << " crosslinks in parallel plane with fibre " <<  (*fibre_neighbor).ID << std::endl;
            if (std::find(this->fibres_crosslinkers.begin(), this->fibres_crosslinkers.end(), (fibre_neighbor)) ==
                this->fibres_crosslinkers.end()) {
                this->fibres_crosslinkers.push_back(fibre_neighbor);
            }
            this->fibres_crosslink_point = point1;
        }
        /*  (3) if fibres are skew (in parallel planes):
            the scalar triple product (P3 - P1) . [(P2 - P1) x (P4 - P3)] != 0
            so intersecting fibres require (P3 - P1) . [(P2 - P1) x (P4 - P3)] = 0
            we include a tolerance on this to allow for fibre radius */
        double test2_tolerance = co_radius;
        double test2 = dot_product(p1_to_p3, FCP);
        if (std::abs(test2) < test2_tolerance) {
            double a = dot_product(p1_to_p2, p1_to_p3) / dot_product(p1_to_p2, p1_to_p2);
            double b = dot_product(p1_to_p2, p3_to_p4) / dot_product(p1_to_p2, p1_to_p2);
            std::vector<double> c(3, 0.0);
            std::vector<double> n(3, 0.0);
            for (int i = 0; i < 3; i++) {
                c[i] = b * p1_to_p2[i] - p3_to_p4[i];
                n[i] = p1_to_p3[i] - a * p1_to_p2[i];
            }
            double t_neighbor = dot_product(c, n) / dot_product(c, c);
            double t_this = a + b * t_neighbor;
            std::vector<double> crosslink_point(3, 0.0);
            for (int i = 0; i < 2; i++) {
                crosslink_point[i] = point1[i] + t_this * p1_to_p2[i];
            }
            /*  For fibres to intersect the "t" values for both line equations
                must lie in [0,1] we include a tolerance to allow for fibre normalized co_radius */
            double tolerance = co_radius/co_length; //(*fibre_neighbor).custom_data["mRadius"] / (2 * fibre->custom_data["mLength"]);
            double lower_bound = 0.0 - tolerance;
            double upper_bound = 1.0 + tolerance;
            if (lower_bound <= t_neighbor && t_neighbor <= upper_bound &&
                lower_bound <= t_this && t_this <= upper_bound) {
                if (std::find(this->fibres_crosslinkers.begin(), this->fibres_crosslinkers.end(), (fibre_neighbor)) ==
                    this->fibres_crosslinkers.end()) {
                    //std::cout << "fibre " << fibre->ID << " crosslinks with fibre " << (*fibre_neighbor).ID << std::endl;
                    this->fibres_crosslinkers.push_back(fibre_neighbor);
                }
                this->fibres_crosslink_point = crosslink_point;
            }
        }
    } else { return; }
}


void PhysiMeSS_Fibre::add_crosslinks() 
{
    for (auto* neighbor : physimess_neighbors)
    {
        if (isFibre(neighbor)) {
            this->check_fibre_crosslinks(static_cast<PhysiMeSS_Fibre*>(neighbor));
        }
    }
}