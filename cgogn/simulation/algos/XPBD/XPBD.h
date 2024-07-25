#ifndef CGOGN_SIMULATION_XPBD_XPBD_H_
#define CGOGN_SIMULATION_XPBD_XPBD_H_

#include <cgogn/core/types/maps/map_base.h>

#include <algorithm>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/algos/volume.h>
#include <cgogn/geometry/types/vector_traits.h>




#define LAME_MU (YOUNG_MODULUS / (2 * (1 + POISSON_RATIO)))
#define LAME_LAMBDA ((YOUNG_MODULUS * POISSON_RATIO) / ((1 + POISSON_RATIO) * (1 - 2 * POISSON_RATIO)))

#define SHEAR_MODULUS (YOUNG_MODULUS / (2 * (1 + POISSON_RATIO)))
#define BULK_MODULUS (YOUNG_MODULUS / (3 * (1 - 2 * POISSON_RATIO)))

#define NUM_SUBSTEP 1

int DENSITY; //Default value
float POISSON_RATIO;
int YOUNG_MODULUS;
float TIME_STEP;


Eigen::Vector3f position_plan;
Eigen::Vector3f dimension;
int angle;
float friction_coefficient;

Eigen::Vector3f pos_cube;
Eigen::Vector3f size_cube;

Eigen::Vector3f bbmin;
Eigen::Vector3f bbmax;


void setDensity (int density_){
    DENSITY=density_;
}

void setPoissonRatio (float poisson_ratio_){
    POISSON_RATIO=poisson_ratio_;
}

void setYoungModule(int young_modulus_){
    YOUNG_MODULUS = young_modulus_;
}

void setTimeStep(float time_step_){
    TIME_STEP = time_step_;
}

//void setPlanPosition(int x_, int y_, int z_, int width_, int height_, int depth_)
//{
//    pos_cube.x() = x_;
//    pos_cube.y() = y_;
//    pos_cube.z() = z_;
//    size_cube.x() = width_;
//    size_cube.y() = height_;
//    size_cube.z() = depth_;
//}

void setPosition(int position_x, int position_y, int position_z)
{
    position_plan.x() = position_x;
    position_plan.y() = position_y;
    position_plan.z() = position_z;
}

void setDimension(int dimension_x, int dimension_y, int dimension_z)
{
    dimension.x() = dimension_x;
    dimension.y() = dimension_y;
    dimension.z() = dimension_z;
}

void setAngle(int angle_){
    angle = angle_;
}
void setFrictionCoef(float friction_coef_){
    friction_coefficient = friction_coef_;
}


std::string csv_data_energy;
std::string csv_data_position;

void Setcsv_energy (std::string filename)
{
    csv_data_energy = filename;
}

void Setcsv_position (std::string filename)
{
    csv_data_position = filename;
}


#define EPS 1e-12

namespace cgogn
{
namespace simulation
{
template <typename MAP>
class XPBD
{
    using Self = XPBD;
    template <typename T>
    using Attribute = typename mesh_traits<MAP>::template Attribute<T>;
    using Vec3 = geometry::Vec3;
    using Mat3d = geometry::Mat3d;
    using Vertex = typename mesh_traits<MAP>::Vertex;
    using Volume = typename mesh_traits<MAP>::Volume;

public:
    // Stable Values
    std::shared_ptr<Attribute<Vec3>> init_pos_;
    std::shared_ptr<Attribute<Vec3>> init_cm_;
    std::shared_ptr<Attribute<double>> init_volume_;
    std::shared_ptr<Attribute<double>> masse_;
    std::shared_ptr<Attribute<Mat3d>> inv_Q_;
    std::shared_ptr<Attribute<bool>> fixed_vertex;
    std::shared_ptr<Attribute<std::vector<Vertex>>> inc_vertices_;
    // Integration Values

    std::shared_ptr<Attribute<Vec3>> pos_;
    std::shared_ptr<Attribute<Vec3>> pos_prev_;
    std::shared_ptr<Attribute<Vec3>> speed_;
    std::shared_ptr<Attribute<Vec3>> f_ext_;
    std::shared_ptr<Attribute<Vec3>> Grad_C_i_;
    std::shared_ptr<Attribute<Vec3>> Grad_C2_i_;

    XPBD()
        : init_pos_(nullptr), init_cm_(nullptr), masse_(nullptr), inc_vertices_(nullptr), pos_(nullptr),
          pos_prev_(nullptr), speed_(nullptr), f_ext_(nullptr), inv_Q_(nullptr), Grad_C_i_(nullptr)
    {
    }




    void init_solver(MAP& m, std::shared_ptr<Attribute<Vec3>> pos)
    {
        pos_ = pos;
        init_pos_ = get_or_add_attribute<Vec3, Vertex>(m, "XPBD_Init_pos");
        init_cm_ = add_attribute<Vec3, Volume>(m, "XPBD_init_cm");

        masse_ = add_attribute<double, Vertex>(m, "XPBD_masse");

        init_volume_ = add_attribute<double, Volume>(m, "XPBD_init_volume");
        inv_Q_ = add_attribute<Mat3d, Volume>(m, "XPBD_inv_Q");
        inc_vertices_ = add_attribute<std::vector<Vertex>, Volume>(m, "XPBD_inc_vertices_vector_");

        pos_prev_ = add_attribute<Vec3, Vertex>(m, "XPBD_pos_prev_");
        speed_ = add_attribute<Vec3, Vertex>(m, "XPBD_speed_");

        f_ext_ = add_attribute<Vec3, Vertex>(m, "XPBD_f_ext");

        Grad_C_i_ = add_attribute<Vec3, Vertex>(m, "XPBD_Grad_C_i");
        Grad_C2_i_ = add_attribute<Vec3, Vertex>(m, "XPBD_Grad_C2_i");

        parallel_foreach_cell(m, [&](Vertex v) -> bool {
            value<Vec3>(m, init_pos_, v) = value<Vec3>(m, pos_, v);
            value<double>(m, masse_, v) = 0;
            value<Vec3>(m, speed_, v) = Vec3(0, 0, 0);
            value<Vec3>(m, f_ext_, v) = Vec3(0, 0, 0);
            return true;
        });
        geometry::compute_volume(m, pos_.get(), init_volume_.get());
        foreach_cell(m, [&](Volume v) -> bool {
            std::vector<Vertex> vertices;
            foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
                vertices.push_back(w);
                return true;
            });
            double masse = value<double>(m, init_volume_, v) * DENSITY / vertices.size();
            for (auto w : vertices)
            {
                value<double>(m, masse_, w) += masse;
            }
            return true;
        });

        parallel_foreach_cell(m, [&](Volume v) -> bool {
            double masse = 0;
            Vec3 cm = Vec3(0, 0, 0);
            foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
                masse += value<double>(m, masse_, w);
                cm += value<double>(m, masse_, w) * value<Vec3>(m, init_pos_, w);
                return true;
            });
            value<Vec3>(m, init_cm_, v) = cm / masse;
            return true;
        });

        parallel_foreach_cell(m, [&](Volume v) -> bool {
            Mat3d Q = Mat3d::Zero();
            foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
                Vec3 init_r_i = value<Vec3>(m, init_pos_, w) - value<Vec3>(m, init_cm_, v);
                double masse = value<double>(m, masse_, w);
                Q += masse * init_r_i * init_r_i.transpose();
                return true;
            });
            value<Mat3d>(m, inv_Q_, v) = Q.inverse();
            return true;
        });
    }

    double energy_calculation(MAP& m, Volume v, double h)
    {
        // Compute center of mass
        double masse_vol = 0;
        Vec3 cm = Vec3(0, 0, 0);
        std::vector<Vertex>& inc_vertices = value<std::vector<Vertex>>(m, inc_vertices_.get(), v);
        for (Vertex w : inc_vertices)
        {
            masse_vol += value<double>(m, masse_, w);
            cm += value<double>(m, masse_, w) * value<Vec3>(m, pos_.get(), w);
        }
        cm /= masse_vol;


        // Compute P
        Mat3d P = Mat3d::Zero();
        for (Vertex w : inc_vertices)
        {
            double m_i = value<double>(m, masse_, w);
            Vec3 r_i = value<Vec3>(m, pos_, w) - cm;
            Vec3 init_r_i = value<Vec3>(m, init_pos_, w) - value<Vec3>(m, init_cm_, v);
            P += m_i * r_i * init_r_i.transpose();
        }

        Mat3d inv_Q = value<Mat3d>(m, inv_Q_, v);


        // compute F = P*Q^-1
        Mat3d F = P * inv_Q;

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                if (i == j)
                {
                    if (fabs(F(i, j) - 1) < EPS)
                        F(i, j) = 1;
                }
                else
                {
                    if (fabs(F(i, j)) < EPS)
                        F(i, j) = 0;
                }
            }
        }

        // Compute PSI_H
        double PSI_H = (LAME_LAMBDA/2) * (F.determinant() - (1 + (LAME_MU / LAME_LAMBDA))) * (F.determinant() - (1 + (LAME_MU / LAME_LAMBDA)));

        // Compute PSI_D
        double PSI_D = (LAME_MU/2) * ((F.transpose() * F).trace() - 3 );


        double PSI_NEO = PSI_H + PSI_D;

        return PSI_NEO;
    }

    Vec3 centre_masse(MAP& m, Volume v, double h)
    {
        double masse_vol = 0;
        Vec3 cm = Vec3(0, 0, 0);
        std::vector<Vertex>& inc_vertices = value<std::vector<Vertex>>(m, inc_vertices_.get(), v);
        for (Vertex w : inc_vertices)
        {
            masse_vol += value<double>(m, masse_, w);
            cm += value<double>(m, masse_, w) * value<Vec3>(m, pos_.get(), w);
        }
        cm /= masse_vol;

        return cm;
    }


    void constraint_Neo_Hookean_H(MAP& m, Volume v, double h)
    {
        // Compute center of mass
        double masse_vol = 0;
        Vec3 cm = Vec3(0, 0, 0);
        std::vector<Vertex>& inc_vertices = value<std::vector<Vertex>>(m, inc_vertices_.get(), v);
        for (Vertex w : inc_vertices)
        {
            masse_vol += value<double>(m, masse_, w);
            cm += value<double>(m, masse_, w) * value<Vec3>(m, pos_.get(), w);
        }
        cm /= masse_vol;

        // Compute P
        Mat3d P = Mat3d::Zero();
        for (Vertex w : inc_vertices)
        {
            double m_i = value<double>(m, masse_, w);
            Vec3 r_i = value<Vec3>(m, pos_, w) - cm;
            Vec3 init_r_i = value<Vec3>(m, init_pos_, w) - value<Vec3>(m, init_cm_, v);
            P += m_i * r_i * init_r_i.transpose();
        }

        Mat3d inv_Q = value<Mat3d>(m, inv_Q_, v);
        // compute F = P*Q^-1
        Mat3d F = P * inv_Q;

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                if (i == j)
                {
                    if (fabs(F(i, j) - 1) < EPS)
                        F(i, j) = 1;
                }
                else
                {
                    if (fabs(F(i, j)) < EPS)
                        F(i, j) = 0;
                }
            }
        }



        // Compute C  = det(F) - (1+MU/LAMBDA)
        double det_F = F.determinant();
        double C = det_F - (1 + LAME_MU / LAME_LAMBDA);

        // Compute Volume
        // double Ve = geometry::volume(m, v, pos_.get());
        double Ve = fabs(det_F * value<double>(m, init_volume_, v));

        // Compute alpha_H = 1/(LAMBDA*V)
        double alpha_h = 1.0f / (LAME_LAMBDA * Ve);

        // Compute denum
        double denum = 0;
        Mat3d T_inv_Q = inv_Q.transpose();
        Mat3d tmp;
        tmp.col(0) = F.col(1).cross(F.col(2));
        tmp.col(1) = F.col(2).cross(F.col(0));
        tmp.col(2) = F.col(0).cross(F.col(1));
        for (Vertex w : inc_vertices)
        {
            // Compute Grad_C_i
            double m_i = value<double>(m, masse_, w);
            Vec3 init_r_i = value<Vec3>(m, init_pos_.get(), w) - value<Vec3>(m, init_cm_, v);
            Vec3 GC = m_i * tmp * T_inv_Q * init_r_i;
            value<Vec3>(m, Grad_C_i_, w) = GC;
            denum += 1.0f / m_i * GC.squaredNorm();
        }
        denum += alpha_h / (h * h);
        // Compute lambda
        double lambda = -C / denum;
        for (Vertex w : inc_vertices)
        {
            if (this->fixed_vertex && value<bool>(m, this->fixed_vertex.get(), w))
            {
                continue;
            }
            double m_i = value<double>(m, masse_, w);
            Vec3 delta_x = lambda * (1 / m_i) * value<Vec3>(m, Grad_C_i_, w);
            value<Vec3>(m, pos_.get(), w) += delta_x;
        }

    }
    void constraint_Neo_Hookean_D(MAP& m, Volume v, double h)
    {
        // Compute center of mass
        double masse_vol = 0;
        Vec3 cm = Vec3(0, 0, 0);
        std::vector<Vertex>& inc_vertices = value<std::vector<Vertex>>(m, inc_vertices_.get(), v);
        for (Vertex w : inc_vertices)
        {
            masse_vol += value<double>(m, masse_, w);
            cm += value<double>(m, masse_, w) * value<Vec3>(m, pos_.get(), w);
        }
        cm /= masse_vol;

        // Compute P
        Mat3d P = Mat3d::Zero();
        for (Vertex w : inc_vertices)
        {
            double m_i = value<double>(m, masse_, w);
            Vec3 r_i = value<Vec3>(m, pos_.get(), w) - cm;
            Vec3 init_r_i = value<Vec3>(m, init_pos_.get(), w) - value<Vec3>(m, init_cm_, v);
            P += m_i * r_i * init_r_i.transpose();
        }

        Mat3d inv_Q = value<Mat3d>(m, inv_Q_, v);
        // compute F = P*Q^-1
        Mat3d F = P * inv_Q;

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                if (i == j)
                {
                    if (fabs(F(i, j) - 1) < EPS)
                        F(i, j) = 1;
                }
                else
                {
                    if (fabs(F(i, j)) < EPS)
                        F(i, j) = 0;
                }
            }
        }


        // Compute C  = sqrt(tr(F^T*F))
        double C = sqrt((F.transpose() * F).trace());

        // Compute Volume
        // double Ve = geometry::volume(m, v, pos_.get());
        double Ve = fabs(F.determinant() * value<double>(m, init_volume_, v));

        // Compute alpha_D = 1/(MU*V)
        double alpha_d = 1.0f / (LAME_MU * Ve);

        // Compute denum
        double denum = 0;
        Mat3d T_inv_Q = inv_Q.transpose();
        double r = sqrt(F.col(0).squaredNorm() + F.col(1).squaredNorm() + F.col(2).squaredNorm());
        for (Vertex w : inc_vertices)
        {
            // Compute Grad_C_i
            double m_i = value<double>(m, masse_, w);
            Vec3 init_r_i = value<Vec3>(m, init_pos_.get(), w) - value<Vec3>(m, init_cm_, v);
            Vec3 GC = m_i / r * F * T_inv_Q * init_r_i;
            value<Vec3>(m, Grad_C_i_, w) = GC;
            denum += 1.0f / m_i * GC.squaredNorm();
        }
        denum += alpha_d / (h * h);
        // Compute lambda
        double lambda = -C / denum;

        for (Vertex w : inc_vertices)
        {
            if (this->fixed_vertex && value<bool>(m, this->fixed_vertex.get(), w))
            {
                continue;
            }
            double m_i = value<double>(m, masse_, w);
            Vec3 delta_x = lambda * (1 / m_i) * value<Vec3>(m, Grad_C_i_, w);
            value<Vec3>(m, pos_.get(), w) += delta_x;
        }

    }
    //*****************Stress tensor***********************************************************//

    double stress_tensor(MAP& m, Volume v)
    {
        // Compute P
        Mat3d P = Mat3d::Zero();
        std::vector<Vertex>& inc_vertices = value<std::vector<Vertex>>(m, inc_vertices_.get(), v);
        Vec3 cm = Vec3(0, 0, 0);
        for (Vertex w : inc_vertices)
        {
            double m_i = value<double>(m, masse_, w);
            Vec3 r_i = value<Vec3>(m, pos_, w) - cm;
            Vec3 init_r_i = value<Vec3>(m, init_pos_, w) - value<Vec3>(m, init_cm_, v);
            P += m_i * r_i * init_r_i.transpose();
        }

        Mat3d inv_Q = value<Mat3d>(m, inv_Q_, v);
        // Compute F
        Mat3d F = P * inv_Q;

        // Compute I

        Mat3d I_ = Mat3d::Zero();
        for(int i=0; i<3; ++i)
            I_(i, i) = 1;

        // Compute ε
        Mat3d strain_tensor = Mat3d::Zero();
        strain_tensor = 0.5*((F.transpose())*F - I_);

        // Compute D*ε
        Mat3d mat_val = Mat3d::Zero();
        mat_val = 2*LAME_MU*strain_tensor + LAME_LAMBDA*strain_tensor.trace()*I_;

        // Compute stress tensor
        Mat3d s_t = Mat3d::Zero();
        s_t = F*mat_val;

//        std::cout << "Stress tensor:" << std::endl;
//        for (int i = 0; i < 3; ++i) {
//            for (int j = 0; j < 3; ++j) {
//                std::cout << s_t(i, j) << " ";
//            }
//            std::cout << std::endl;
//        }

        double sigma_vm = (1/sqrt(2)) * sqrt(  (s_t(1, 1)-s_t(2, 2)) * (s_t(1, 1)-s_t(2, 2))
                                           +   (s_t(2, 2)-s_t(3, 3)) * (s_t(2, 2)-s_t(3, 3))
                                           +   (s_t(3, 3)-s_t(1, 1)) * (s_t(3, 3)-s_t(1, 1))
                                           +   6 * (s_t(1, 2)*s_t(1, 2) + s_t(1, 3)*s_t(1, 3) + s_t(2, 3)*s_t(2, 3)
                                           )
                                           );

//        std::cout << "Von Mises = " << sigma_vm << std::endl;

        return sigma_vm;
    }


    void constraint_Zero_Energy(MAP& m, Volume v, double)
    {
        // Compute center of mass
        double masse_vol = 0;
        Vec3 cm = Vec3(0, 0, 0);
        std::vector<Vertex>& inc_vertices = value<std::vector<Vertex>>(m, inc_vertices_.get(), v);
        for (Vertex w : inc_vertices)
        {
            masse_vol += value<double>(m, masse_, w);
            cm += value<double>(m, masse_, w) * value<Vec3>(m, pos_.get(), w);
        }
        cm /= masse_vol;

        // Compute P
        Mat3d P = Mat3d::Zero();
        for (Vertex w : inc_vertices)
        {
            double m_i = value<double>(m, masse_, w);
            Vec3 r_i = value<Vec3>(m, pos_.get(), w) - cm;
            Vec3 init_r_i = value<Vec3>(m, init_pos_.get(), w) - value<Vec3>(m, init_cm_, v);
            P += m_i * r_i * init_r_i.transpose();
        }

        Mat3d inv_Q = value<Mat3d>(m, inv_Q_, v);
        // compute F = P*Q^-1
        Mat3d F = P * inv_Q;

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                if (i == j)
                {
                    if (fabs(F(i, j) - 1) < EPS)
                        F(i, j) = 1;
                }
                else
                {
                    if (fabs(F(i, j)) < EPS)
                        F(i, j) = 0;
                }
            }
        }

        for (Vertex w : inc_vertices)
        {
            if (this->fixed_vertex && value<bool>(m, this->fixed_vertex.get(), w))
            {
                continue;
            }
            Vec3 init_r_i = value<Vec3>(m, init_pos_.get(), w) - value<Vec3>(m, init_cm_, v);
            value<Vec3>(m, pos_.get(), w) = cm + F * init_r_i;
        }
    }
    void applyDamping(MAP& m, Volume v, double damping_coeff, double time_step)
    {
        Vec3 x_cm = Vec3::Zero();
        Vec3 v_cm = Vec3::Zero();
        Vec3 L = Vec3::Zero();
        Mat3d I = Mat3d::Zero();
        double sm_i = 0;
        foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
            double m_i = value<double>(m, masse_, w);
            x_cm += m_i * value<Vec3>(m, pos_.get(), w);
            v_cm += m_i * value<Vec3>(m, speed_, w);
            sm_i += m_i;
            return true;
        });
        x_cm /= sm_i;
        v_cm /= sm_i;
        foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
            Vec3 r_i = value<Vec3>(m, pos_.get(), w) - x_cm;
            double m_i = value<double>(m, masse_, w);
            L += r_i.cross(m_i * value<Vec3>(m, speed_, w));
            Mat3d Ri;
            Ri << 0, -r_i(2), r_i(1), r_i(2), 0, -r_i(0), -r_i(1), r_i(0), 0;
            I += Ri * Ri.transpose() * m_i;
            return true;
        });
        Vec3 omega = I.inverse() * L;
        foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
            Vec3 r_i = value<Vec3>(m, pos_.get(), w) - x_cm;
            Vec3 new_v_i = v_cm + omega.cross(r_i);
            value<Vec3>(m, speed_, w) +=
                std::min(damping_coeff * time_step, 1.0) * (new_v_i - value<Vec3>(m, speed_, w));
            return true;
        });
    }


#define SHOW_PERFORMANCE_LOG 1
    void solver(MAP& m, double timestep)
    {
        std::clock_t start;
        double duration;
        start = std::clock();
        double h = timestep / NUM_SUBSTEP;
        std::vector<Volume> vec_volume;
        std::vector<Vertex> vec_vertices;
        foreach_cell(m, [&](Vertex v) -> bool {
            vec_vertices.push_back(v);
            return true;
        });
        foreach_cell(m, [&](Volume v) -> bool {
            vec_volume.push_back(v);
            value<std::vector<Vertex>>(m, inc_vertices_.get(), v).clear();
            foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
                value<std::vector<Vertex>>(m, inc_vertices_.get(), v).push_back(w);
                return true;
            });
            return true;
        });
        duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
#if SHOW_PERFORMANCE_LOG
        std::cout << "\033[1;32mtime init XPBD : \033[0m" << duration << std::endl;
#endif
        std::srand(455279423);
        std::random_shuffle(vec_volume.begin(),vec_volume.end());
        start = std::clock();
        for (int i = 0; i < NUM_SUBSTEP; i++)
        {
            // Initialisation sub step
            for (Vertex v : vec_vertices)
            {
                if (this->fixed_vertex && value<bool>(m, this->fixed_vertex.get(), v))
                {
                    continue;
                }
                value<Vec3>(m, pos_prev_, v) = value<Vec3>(m, pos_, v);
                value<Vec3>(m, speed_, v) += h * value<Vec3>(m, f_ext_, v) / value<double>(m, masse_, v);
                value<Vec3>(m, pos_, v) += h * value<Vec3>(m, speed_, v);
            }
            // Constraint solver
            for (Volume v : vec_volume)
            {
                constraint_Neo_Hookean_H(m, v, h);
                constraint_Neo_Hookean_D(m, v, h);
                constraint_Zero_Energy(m, v, h);
            }


            // Speed Solver
            for (Vertex v : vec_vertices)
            {
                if (this->fixed_vertex && value<bool>(m, this->fixed_vertex.get(), v))
                {
                    continue;
                }
                Vec3 new_v = (value<Vec3>(m, pos_, v) - value<Vec3>(m, pos_prev_, v)) / h;
                for (int i = 0; i < 3; i++)
                    if (fabs(new_v[i]) < EPS)
                        new_v[i] = 0;
                // value<Vec3>(m, speed_, v) = (1 - (0.005 * h)) * new_v;
                value<Vec3>(m, speed_, v) = new_v;
            }
            // Damping
            foreach_cell(m, [&](Volume v) -> bool {
                applyDamping(m, v, 0.1, timestep);
                return true;
            });
        }

        duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
#if SHOW_PERFORMANCE_LOG
        std::cout << "\033[1;36mtime resolve XPBD : \033[0m" << duration << std::endl;
#endif
    }
};
} // namespace simulation
} // namespace cgogn
#endif // CGOGN_SIMULATION_XPBD_XPBD_H_
