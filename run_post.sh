#!/bin/bash

# Check if a path is provided
if [ -z "$1" ]; then
    echo "Error: You must provide the path to the OpenFOAM case."
    echo "Usage: ./run_post.sh /path/to/your/case"
    exit 1
fi

# Convert to absolute path and export it for Python
export CASE_DIR=$(realpath "$1")

# Check if the case directory exists
if [ ! -d "$CASE_DIR" ]; then
    echo "Error: Case folder $CASE_DIR does not exist."
    exit 1
fi

echo "Processing case: $CASE_DIR"
cd "$CASE_DIR" || exit

# Handle parallel reconstruction (only latest time)
if ls -d processor* >/dev/null 2>&1; then
    echo "Parallel processor folders detected. Reconstructing latest time step..."
    reconstructPar -latestTime
    
    # Check if reconstruction was successful before deleting
    if [ $? -eq 0 ]; then
        echo "Reconstruction successful. Deleting processor folders..."
        rm -rf processor*
    else
        echo "Error during reconstructPar. Stopping."
        exit 1
    fi
else
    echo "No processor folders detected. Skipping reconstruction."
fi

# Handle VTK generation (only latest time)
if [ ! -d "VTK" ]; then
    echo "VTK folder not found. Running foamToVTK for the latest time step..."
    foamToVTK -latestTime
else
    echo "VTK folder already exists. Skipping foamToVTK."
fi

# Return to the original directory where the script was called
cd - > /dev/null
TOOL_DIR=$(pwd)

CASE_NAME=$(basename "$CASE_DIR")
PROCESSED_DIR="$TOOL_DIR/Processed_cases/${CASE_NAME}_processed"

mkdir -p "$PROCESSED_DIR/raw_data"
mkdir -p "$PROCESSED_DIR/notebooks"
mkdir -p "$PROCESSED_DIR/figures"

export RAW_DATA_DIR="$PROCESSED_DIR/raw_data"
export FIGURES_DIR="$PROCESSED_DIR/figures"


NETWORK_FILE="$RAW_DATA_DIR/junction_links.txt"
COORDS_FILE="$RAW_DATA_DIR/junction_coordinates.txt"

if [ -f "$NETWORK_FILE" ]; then
    echo "FAST MODE: Network files already detected!"
    echo "Skipping the 7-minute Jupyter notebook execution."
else
    echo "Activating 'tubes' virtual environment..."
    source "$TOOL_DIR/tubes/bin/activate"

    echo "Running Jupyter notebook automatically..."
    jupyter nbconvert --to notebook --execute --output-dir="$PROCESSED_DIR/notebooks" --output="${CASE_NAME}_halfwidths.ipynb" halfwidths_and_fluxes_v6.ipynb
    

    if [ $? -eq 0 ]; then
        echo "Post-processing complete! Files saved in $PROCESSED_DIR ."
        deactivate
    else
        echo "ERROR: Jupyter Notebook execution failed. Check the logs above."
        deactivate
        exit 1
    fi
fi


echo "Running C++ scripts..."

CPP_SOURCE="random_splitting_def.cpp myfunctions.cpp"
EXEC_NAME="./random_splitting"

# compile
echo "Compiling $CPP_SOURCE..."
g++ -std=c++17 -O3 -I"$TOOL_DIR/eigen-5.0.0" $CPP_SOURCE -o "$EXEC_NAME"

if [ $? -ne 0 ]; then
    echo "FATAL ERROR: C++ compilation failed."
    exit 1
fi

# Execute
N_REALS=1000   
N_BINS=50      

echo "Executing..."
$EXEC_NAME "$NETWORK_FILE" $N_REALS $N_BINS "$COORDS_FILE"

if [ $? -eq 0 ]; then
    if [ -d "output_def" ]; then
        mv output_def/* "$RAW_DATA_DIR/"
        rmdir output_def
    fi
    echo "SUCCESS: C++ Analysis complete! All files organized in:"
    echo "$PROCESSED_DIR"

    echo "Activating 'tubes' virtual environment for visualization..."
    source "$TOOL_DIR/tubes/bin/activate"

    echo "Running first visualization notebook..."
    jupyter nbconvert --to notebook --execute --output-dir="$PROCESSED_DIR/notebooks" --output="${CASE_NAME}_visualization.ipynb" random_splitting_well_plotted_def.ipynb

    if [ $? -eq 0 ]; then
        echo "Success."
    else
        echo "Visualization notebook execution failed"
        deactivate
        exit 1
    fi
    deactivate

    echo "Running second C++ script (Pore Tests)..."

    CPP_SOURCE_2="pore_tests.cpp myfunctions.cpp"
    EXEC_NAME_2="./pore_tests_exec"

    # compile
    echo "Compiling $CPP_SOURCE_2..."
    g++ -std=c++17 -O3 -I"$TOOL_DIR/eigen-5.0.0" $CPP_SOURCE_2 -o "$EXEC_NAME_2"

    if [ $? -ne 0 ]; then
        echo "FATAL ERROR: Second C++ compilation (pore_tests) failed."
        exit 1
    fi

    # Execute
    echo "Executing pore_tests..."
    $EXEC_NAME_2 "$RAW_DATA_DIR" "$NETWORK_FILE" "$COORDS_FILE"

    if [ $? -eq 0 ]; then
        echo "SUCCESS: Pore tests complete! All .dat histograms saved to raw_data."
    else
        echo "FATAL ERROR: pore_tests execution failed."
        exit 1
    fi

    echo "Running third C++ script (Variance Scaling)..."

    CPP_SOURCE_3="variance_scaling.cpp myfunctions.cpp"
    EXEC_NAME_3="./variance_scaling_exec"

    # compile
    echo "Compiling $CPP_SOURCE_3..."
    g++ -std=c++17 -O3 -I"$TOOL_DIR/eigen-5.0.0" $CPP_SOURCE_3 -o "$EXEC_NAME_3"

    if [ $? -ne 0 ]; then
        echo "Third C++ compilation (variance_scaling) failed."
        exit 1
    fi

    # Execute
    echo "Executing variance_scaling..."
    $EXEC_NAME_3 "$RAW_DATA_DIR" "$NETWORK_FILE" "$COORDS_FILE"

    if [ $? -eq 0 ]; then
        echo "SUCCESS: Variance scaling complete! All scaling .dat files saved to raw_data."
    else
        echo "FATAL ERROR: variance_scaling execution failed."
        exit 1
    fi

    
    echo "Activating 'tubes' virtual environment for final visualization..."
    source "$TOOL_DIR/tubes/bin/activate"

    echo "Running variance scaling visualization notebook..."
    jupyter nbconvert --to notebook --execute \
        --output-dir="$PROCESSED_DIR/notebooks" \
        --output="${CASE_NAME}_variance_scaling_plotted.ipynb" \
        variance_scaling_well_plotted.ipynb

    if [ $? -eq 0 ]; then
        echo "SUCCESS: Final visualization complete! PDF saved in figures."
        echo "🏆 THE GRAND FINALE: ALL PROCESSES FINISHED SUCCESSFULLY!"
       
    else
        echo "ERROR: Final visualization notebook execution failed."
        deactivate
        exit 1
    fi
    deactivate
    
    


else
    echo "FATAL ERROR: C++ execution failed."
    exit 1
fi