#!/usr/bin/env julia

"""
Simple script to compare new algorithm output with reference data.
Usage: julia scripts/compare_with_reference.jl
"""

using JSON

function compare_with_reference(models, ref_file="14bus_ref_sol.json"; tolerance=1e-5)
    """
    Compare current bus values with reference values.
    
    Args:
        models: Current models dictionary with bus data
        ref_file: Reference JSON file path
        tolerance: Tolerance for comparison
    
    Returns:
        true if all values match within tolerance, false otherwise
    """
    
    # Check if reference file exists
    if !isfile(ref_file)
        println("Reference file '$ref_file' not found!")
        println("Make sure to run the power flow simulation first to generate the reference.")
        return false
    end
    
    # Load reference data
    ref_data = JSON.parsefile(ref_file)
    
    # Convert JSON arrays to Float64 arrays
    v_ref = Float64.(ref_data["v"])
    theta_ref = Float64.(ref_data["theta"])
    vd_ref = Float64.(ref_data["vd"])
    vq_ref = Float64.(ref_data["vq"])
    
    # Compare values
    v_match = all(abs.(models.bus.v - v_ref) .< tolerance)
    theta_match = all(abs.(models.bus.theta - theta_ref) .< tolerance)
    vd_match = all(abs.(models.bus.vd - vd_ref) .< tolerance)
    vq_match = all(abs.(models.bus.vq - vq_ref) .< tolerance)
    
    # Print results
    println("=== Comparison Results ===")
    println("Reference file: $ref_file")
    println("Tolerance: $tolerance")
    println()
    println("v match:     $(v_match ? "✓" : "✗")")
    println("theta match: $(theta_match ? "✓" : "✗")")
    println("vd match:    $(vd_match ? "✓" : "✗")")
    println("vq match:    $(vq_match ? "✓" : "✗")")
    println()
    
    if v_match && theta_match && vd_match && vq_match
        println("✅ All values match reference within tolerance!")
        return true
    else
        println("❌ Some values differ from reference!")
        
        # Show differences if any
        if !v_match
            println("\nv differences:")
            for i in eachindex(models.bus.v)
                diff = abs(models.bus.v[i] - v_ref[i])
                if diff >= tolerance
                    println("  Bus $(models.bus.idx[i]): $(models.bus.v[i]) vs $(v_ref[i]) (diff: $diff)")
                end
            end
        end
        
        if !theta_match
            println("\ntheta differences:")
            for i in eachindex(models.bus.theta)
                diff = abs(models.bus.theta[i] - theta_ref[i])
                if diff >= tolerance
                    println("  Bus $(models.bus.idx[i]): $(models.bus.theta[i]) vs $(theta_ref[i]) (diff: $diff)")
                end
            end
        end
        
        if !vd_match
            println("\nvd differences:")
            for i in eachindex(models.bus.vd)
                diff = abs(models.bus.vd[i] - vd_ref[i])
                if diff >= tolerance
                    println("  Bus $(models.bus.idx[i]): $(models.bus.vd[i]) vs $(vd_ref[i]) (diff: $diff)")
                end
            end
        end
        
        if !vq_match
            println("\nvq differences:")
            for i in eachindex(models.bus.vq)
                diff = abs(models.bus.vq[i] - vq_ref[i])
                if diff >= tolerance
                    println("  Bus $(models.bus.idx[i]): $(models.bus.vq[i]) vs $(vq_ref[i]) (diff: $diff)")
                end
            end
        end
        
        return false
    end
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    println("=== Bus Values Comparison Script ===")
    println()
    
    # Check if reference file exists
    ref_file = "14bus_ref_sol.json"
    if !isfile(ref_file)
        println("❌ Reference file '$ref_file' not found!")
        println()
        println("To use this script:")
        println("1. First run your power flow simulation to generate the reference")
        println("2. Then run this script to compare new results")
        println()
        println("Example:")
        println("  julia src/power_flow_traditional.jl  # Generate reference")
        println("  julia scripts/compare_with_reference.jl  # Compare results")
        exit(1)
    end
    
    println("Reference file found: $ref_file")
    println("Note: This script expects 'models' to be available in the current scope.")
    println("Make sure to run this after your power flow simulation.")
    println()
    
    # Try to load reference data to show what's expected
    ref_data = JSON.parsefile(ref_file)
    println("Reference data contains:")
    println("  Number of buses: $(length(ref_data["bus_indices"]))")
    println("  Bus indices: $(ref_data["bus_indices"])")
    println()
    println("To compare your results, call:")
    println("  compare_with_reference(models)")
    println()
    println("Or with custom tolerance:")
    println("  compare_with_reference(models, tolerance=1e-8)")
end
