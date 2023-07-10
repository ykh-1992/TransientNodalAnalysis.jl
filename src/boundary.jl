struct IARF <: Boundary end

IARF(_) = IARF()

dimensionless_well_pressure(::IARF, t_D, s, _) = -0.5*ei(-0.25/t_D) + s

struct NoFlow <: Boundary
    distance::Float64 # ft, distance to boundary
end

function dimensionless_well_pressure(boundary::NoFlow, t_D, s, r_w)
    r_eD = boundary.distance/r_w # ft, dimensionless radius to reservoir extent
    p_wD = -0.5*ei(-0.25/t_D) - 0.5*ei(-r_eD^2/t_D) + s # -, dimensionless well pressure
    return p_wD
end

struct ParallelNoFlow <: Boundary
    distance::Float64 # ft, distance to boundary
end

function dimensionless_well_pressure(boundary::ParallelNoFlow, t_D, s, r_w)
    r_eD = boundary.distance/r_w      # ft, dimensionless radius to reservoir extent
    N_iw = num_image_wells(t_D, r_eD) # -, number of image wells on one side
    
    p_wD = 0.0 # -, dimensionless well pressure
    for i in 1:N_iw
        p_wD += ei(-(i*r_eD)^2/t_D) # -, dimensionless well pressure
    end

    p_wD = -0.5*ei(-0.25/t_D) - p_wD + s # -, dimensionless well pressure

    return p_wD
end

struct BoxNoFlow <: Boundary
    distance::Float64 # ft, distance to boundary
end

function dimensionless_well_pressure(boundary::BoxNoFlow, t_D, s, r_w)
    r_eD = boundary.distance/r_w      # ft, dimensionless radius to reservoir extent
    N_iw = num_image_wells(t_D, r_eD) # -, number of image wells on one side

    p_wD_axis_diag = 0.0 # -, dimensionless well pressure for axis and diagonal
    for i in 1:N_iw
        x = -(i*r_eD)^2/t_D
        p_wD_axis_diag += ei(x) + ei(2*x) # -
    end

    p_wD_between = 0.0 # -, dimensionless well pressure for in-between
    for i in 2:N_iw, j in 1:(i - 1)
        p_wD_between += ei(-(i^2 + j^2)*r_eD^2/t_D) # -
    end

    p_wD = -0.5*ei(-0.25/t_D) - 2*p_wD_axis_diag - 4*p_wD_between + s # -, dimensionless well pressure

    return p_wD
end

struct ConstantPressure <: Boundary
    distance::Float64 # ft, distance to boundary
end

function dimensionless_well_pressure(boundary::ConstantPressure, t_D, s, r_w)
    r_eD = boundary.distance/r_w # ft, dimensionless radius to reservoir extent
    p_wD = -0.5*ei(-0.25/t_D) + 0.5*ei(-r_eD^2/t_D) + s # -, dimensionless well pressure
    return p_wD
end

struct ParallelConstantPressure <: Boundary
    distance::Float64 # ft, distance to boundary
end

function dimensionless_well_pressure(boundary::ParallelConstantPressure, t_D, s, r_w)
    r_eD = boundary.distance/r_w      # ft, dimensionless radius to reservoir extent
    N_iw = num_image_wells(t_D, r_eD) # -, number of image wells on one side
    
    p_wD = 0.0 # -, dimensionless well pressure
    for i in 1:N_iw
        p_wD += (-1)^i*ei(-(i*r_eD)^2/t_D) # -, dimensionless well pressure
    end

    p_wD = -0.5*ei(-0.25/t_D) - p_wD + s # -, dimensionless well pressure

    return p_wD
end

struct BoxConstantPressure <: Boundary
    distance::Float64 # ft, distance to boundary
end

function dimensionless_well_pressure(boundary::BoxConstantPressure, t_D, s, r_w)
    r_eD = boundary.distance/r_w      # ft, dimensionless radius to reservoir extent
    N_iw = num_image_wells(t_D, r_eD) # -, number of image wells on one side

    p_wD_axis_diag = 0.0 # -, dimensionless well pressure for axis and diagonal
    for i in 1:N_iw
        x = -(i*r_eD)^2/t_D
        p_wD_axis_diag += (-1)^i*ei(x) + ei(2*x) # -
    end

    p_wD_between = 0.0 # -, dimensionless well pressure for in-between
    for i in 2:N_iw, j in 1:(i - 1)
        p_wD_between += (-1)^(i + j)*ei(-(i^2 + j^2)*r_eD^2/t_D) # -
    end

    p_wD = -0.5*ei(-0.25/t_D) - 2*p_wD_axis_diag - 4*p_wD_between + s # -, dimensionless well pressure

    return p_wD
end

num_image_wells(t_D, r_eD) = ceil(Int, 2*sqrt(t_D)/r_eD)

ei(x) = expinti(x)
