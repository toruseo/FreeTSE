import os
import pytest # Assuming pytest is or will be the test runner
from FreeTSE import FreeTSE # Import directly from package now
import importlib.resources
import contextlib

# We might need to adjust sys.path if FreeTSE.py is not directly importable
# For now, assume it's in the same directory or PYTHONPATH is set up.

def test_ngsim_grid_scenario():
    output_dir = "test_output"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    tse = FreeTSE(gui_mode=False)
    assert tse is not None, "FreeTSE object should be created."
    assert tse.gui_mode == False, "gui_mode should be False."

    with contextlib.ExitStack() as stack:
        speed_data_path_obj = stack.enter_context(importlib.resources.as_file(importlib.resources.files('FreeTSE.dat').joinpath('ngsim_grid_speed_all.csv')))
        flow_data_path_obj = stack.enter_context(importlib.resources.as_file(importlib.resources.files('FreeTSE.dat').joinpath('ngsim_grid_flow_200m.csv')))
        
        speed_data_name = str(speed_data_path_obj)
        flow_data_name = str(flow_data_path_obj)

        tse.set_scenario(
            name="ngsim_trajectory",
            dt=4,
            dx=100,
            mint=0,
            maxt=800,
            minx=0,
            maxx=500,
            number_of_lanes=5,
            speed_data_name=speed_data_name,
            speed_label_t="t",
            speed_label_x="x",
            speed_label_v="v",
            density_data_name=None,
            density_label_t="t", # Default, not used if density_data_name is None
            density_label_x="x", # Default, not used if density_data_name is None
            density_label_k="k", # Default, not used if density_data_name is None
            flow_data_name=flow_data_name,
            flow_label_t="t",
            flow_label_x="x",
            flow_label_q="q",
            density_dat_true_name=None, # No GroundTruth in ngsim_grid.ini
            true_density_label_t="t", # Default
            true_density_label_x="x", # Default
            true_density_label_k="k", # Default
            flow_dat_true_name=None, # No GroundTruth in ngsim_grid.ini
            true_flow_label_t="t", # Default
            true_flow_label_x="x", # Default
            true_flow_label_q="q"  # Default
        )
        assert tse.name == "ngsim_trajectory", "Scenario name should be set."
        assert tse.dt == 4, "dt should be set."
        assert tse.dx == 100, "dx should be set."
        assert tse.speed_data_name == speed_data_name, "Speed data name should be set."
        assert tse.flow_data_name == flow_data_name, "Flow data name should be set."
        assert tse.groundtruth == False, "Groundtruth should be False for this scenario."
        # Check that data loading part of set_data (called by set_scenario) did something
        assert hasattr(tse, 'vv'), "Speed data (vv) should be loaded."
        assert hasattr(tse, 'kk'), "Density/flow data (kk) should be loaded."
        assert tse.vv is not None, "Speed data (vv) should not be None."
        assert tse.kk is not None, "Density/flow data (kk) should not be None."
        tse.estimation()
        assert hasattr(tse, 'k_smoo'), "Smoothed density (k_smoo) should exist after estimation."
        assert tse.k_smoo is not None, "Smoothed density (k_smoo) should not be None after estimation."
        assert hasattr(tse, 'N'), "Cumulative curves (N) should exist after estimation."
        assert tse.N is not None, "Cumulative curves (N) should not be None after estimation."
        tse.accuracy_evaluation() # Should run fine even with self.groundtruth = False
        # Since groundtruth is False, RMSE/MAPE attributes might not be set or set to default.
        # The main check is that it runs without error.
        # If it sets them to None or 0, we could assert that. For now, no specific assertion here beyond no-error.
        pass # Covered by pytest not catching an exception

        test_output_filename = os.path.join(output_dir, "test_output_ngsim_grid.csv")
        tse.save_results(test_output_filename)

        # For visualize, replicate flags based on plan, explicitly set save=False
        tse.visualize(
            smooth=1,       # tsd_mean
            true=0,         # tsd_true (groundtruth is False)
            speed=1,        # tsd_speed
            timeseries=1,   # timeseries
            inputdata=1,    # tsd_data
            scatter=1,      # scat
            save=True,     # Changed to True
            fname=os.path.join(output_dir, "viz_ngsim_grid") # Changed fname
        )

        # Placeholder for assertions to be added in the next step
        # For now, check if the output file was created
        assert os.path.exists(test_output_filename), f"Output file {test_output_filename} was not created."
