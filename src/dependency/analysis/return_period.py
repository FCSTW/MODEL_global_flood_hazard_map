import numpy as np
import os
import itertools as it
import rasterio
from . import analysis_process
from .. import environment as env

env.init()

class ReturnPeriodAnalysis(analysis_process.AnalysisProcess):

    def __init__(self):

        super().__init__()

        # ================================================================================
        #
        # Update the attributes
        #
        # ================================================================================
        
        # Read the reference raster file and update the attributes
        self.read_reference_raster_file()

        # ================================================================================
        #
        # Call the methods to calculate the exceedance probability
        #
        # ================================================================================

        self.calc_exceedance_probability()

        return
    
    def read_reference_raster_file(self) -> None:

        """
        This method reads the reference raster files and performs the preprocessing procedures.
        ================================================================================

        Arguments:

            None

        Returns:

            None
        """

        # Read the reference raster file
        raster_array, file_crs, file_transform = self._read_raster_data()

        # ================================================================================
        #
        # Update the attributes
        #
        # ================================================================================

        self.ref_shape     = raster_array.shape[1:]
        self.ref_crs       = file_crs
        self.ref_transform = file_transform

        return
    
    def calc_exceedance_probability(self) -> None:

        """
        This method reads the raster files and calculates the exceedance probability.
        ================================================================================

        Arguments:

            None

        Returns:

            None
        """

        # ================================================================================
        #
        # Read the raster files
        #
        # ================================================================================

        # Read inuncoast files
        # Loop over climate scenarios, years, subsidence, and projections

        for climate_scenario, year in it.product(self.list_climate_scenario, self.list_year):
            
            # ================================================================================
            #
            # Skip the loop if the combination of climate scenario and year is not valid
            #
            # ================================================================================

            if ((climate_scenario == 'historical') and (year != 'hist')) or ((climate_scenario != 'historical') and (year == 'hist')):
                
                continue

            for subsidence, projection in it.product(self.list_subsidence, self.list_projection):
                
                # ================================================================================
                #
                # Skip the loop if the combination of climate scenario and project is not valid
                #
                # ================================================================================
                    
                if ((climate_scenario == 'historical') and (projection != '0')):
                    
                    continue

                # ================================================================================
                #
                # Read the raster files and output the results
                #
                # ================================================================================
                
                # Print message
                print(f'Reading the raster files of inuncoast: {climate_scenario}, {year}, {subsidence}, {projection}')

                # Set the output file
                output_file_flooding_rp_idx = f'flooding_rp_idx.inuncoast.{climate_scenario}.{subsidence}.{year}.{projection}.tif'

                # Check if the output file exists
                if (os.path.isfile(os.path.join(env.CONFIG['path']['output'], 'flooding_rp_idx_rawdata', output_file_flooding_rp_idx))):
                    
                    continue
                
                # Create a array to store the raster data
                raster_data_multi_rp = np.full((len(self.list_returnperiod), self.ref_shape[0], self.ref_shape[1]), np.nan, dtype=np.int8)

                # ================================================================================
                #
                # Read the raster files of different return periods
                #
                # ================================================================================

                for return_period in self.list_returnperiod:

                    # Set the file name
                    source_file_name = 'inuncoast_{climatescenario}_{subsidence}_{year}_{returnperiod}_{projection}.tif'.format(
                        climatescenario=climate_scenario,
                        subsidence=subsidence,
                        year=year,
                        returnperiod=return_period,
                        projection=projection,
                    )

                    # Read the raster file
                    raster_array, _, _ = self._read_raster_data(source_file_name)

                    # Check if the file exists
                    if (raster_array is None):

                        continue

                    # Print message
                    print(f'Reading the file {source_file_name}')

                    # Update the raster data
                    raster_data_multi_rp[self.list_returnperiod.index(return_period), ...] = raster_array.squeeze()
                
                # Calculate the dummy of the flooding
                flooding_rp_idx = self._calc_flooding_stats(raster_data_multi_rp)
                
                # ================================================================================
                #
                # Output the raster file
                #
                # ================================================================================

                # Output the raster file
                self._output_raster_file(
                    output_file_flooding_rp_idx,
                    flooding_rp_idx,
                    'flooding_rp_idx',
                    self.ref_crs,
                    self.ref_transform,
                )

        # Read inunriver files
        # Loop over climate scenarios, years and GCMs
        
        for climate_scenario, year in it.product(self.list_climate_scenario, self.list_year):
            
            # ================================================================================
            #
            # Skip the loop if the combination of climate scenario and year is not valid
            #
            # ================================================================================

            if ((climate_scenario == 'historical') and (year != 'hist')) or ((climate_scenario != 'historical') and (year == 'hist')):
                
                continue

            for model in self.list_model:

                # ================================================================================
                #
                # Skip the loop if the combination of climate scenario and model is not valid
                #
                # ================================================================================
                    
                if ((climate_scenario == 'historical') and (model != '000000000WATCH')) or ((climate_scenario != 'historical') and (model == '000000000WATCH')):
                    
                    continue

                # ================================================================================
                #
                # Read the raster files and output the results
                #
                # ================================================================================
                
                # Print message
                print(f'Reading the raster files of inunriver: {climate_scenario}, {year}, {model}')

                # Set the output file
                output_file_flooding_rp_idx = f'flooding_rp_idx.inunriver.{climate_scenario}.{model}.{year}.tif'

                # Check if the output file exists
                if (os.path.isfile(os.path.join(env.CONFIG['path']['output'], 'flooding_rp_idx_rawdata', output_file_flooding_rp_idx))):
                    
                    continue
                
                # Create a array to store the raster data
                raster_data_multi_rp = np.full((len(self.list_returnperiod), self.ref_shape[0], self.ref_shape[1]), np.nan, dtype=np.int8)

                # ================================================================================
                #
                # Read the raster files of different return periods
                #
                # ================================================================================

                for return_period in self.list_returnperiod:

                    # Set the file name
                    source_file_name = 'inunriver_{climatescenario}_{model}_{year}_{returnperiod}.tif'.format(
                        climatescenario=climate_scenario,
                        year=year,
                        returnperiod=return_period,
                        model=model,
                    )

                    # Read the raster file
                    raster_array, _, _ = self._read_raster_data(source_file_name)

                    # Check if the file exists
                    if (raster_array is None):

                        continue

                    # Print message
                    print(f'Reading the file {source_file_name}')

                    # Aggregate the model outputs and update the raster data
                    raster_data_multi_rp[self.list_returnperiod.index(return_period), ...] = raster_array.squeeze()
            
                # Calculate the dummy of the flooding
                flooding_rp_idx = self._calc_flooding_stats(raster_data_multi_rp)
                
                # ================================================================================
                #
                # Output the raster file
                #
                # ================================================================================

                # Output the raster file
                self._output_raster_file(
                    output_file_flooding_rp_idx,
                    flooding_rp_idx,
                    'flooding_rp_idx',
                    self.ref_crs,
                    self.ref_transform,
                )

        return
    
    def _read_raster_data(self, source_file_name: str=None) -> tuple[np.ndarray, rasterio.crs.CRS, rasterio.transform.Affine]:

        """
        This method processes the read raster files.
        ================================================================================

        Arguments:

            source_file_name (str): The name of the source file. If the argument is not provided, the first file in the directory will be read.

        Returns:

            raster_array (np.ndarray): The array of the raster file.

            file_crs (rasterio.crs.CRS): The coordinate reference system of the raster file.

            file_transform (rasterio.transform.Affine): The affine transformation of the raster file.
        """

        # If the file name is not provided, read the first file in the directory
        if (source_file_name is None):

            source_file_name = os.listdir(env.CONFIG['path']['data'])[0]

        # ================================================================================
        #
        # Read the raster file
        #
        # ================================================================================

        # Check if the file exists
        if not (os.path.isfile(os.path.join(env.CONFIG['path']['data'], source_file_name))):

            print(f'The file {source_file_name} does not exist.')
            
            return None, None, None

        # Read the raster file
        with rasterio.open(os.path.join(env.CONFIG['path']['data'], source_file_name)) as f:

            # Get the properties of the raster file
            file_crs=f.crs
            file_transform=f.transform

            # Read the raster file
            raster_array = f.read()

        # Replace the no data value with NaN
        raster_array = np.where(raster_array == -9999, 0, raster_array)

        # Scale the float data to 10 times integer
        raster_array = (raster_array * 10).astype(np.int8)

        return raster_array, file_crs, file_transform
    
    def _output_raster_file(
            self,
            output_file_name: str,
            raster_array: np.ndarray,
            variable_name: str,
            file_crs: rasterio.crs.CRS,
            file_transform: rasterio.transform.Affine
        ) -> None:

        """
        This method processes the output raster files.
        ================================================================================

        Arguments:

            output_file_name (str): The name of the output file.

            raster_array (np.ndarray): The array of the raster file.

            variable_name (str): The name of the variable.

            file_crs (rasterio.crs.CRS): The coordinate reference system of the raster file.

            file_transform (rasterio.transform.Affine): The affine transformation of the raster file.

        Returns:

            None
        """

        # Print message
        print(f'Outputting the file {output_file_name}')

        output_path = os.path.join(env.CONFIG['path']['output'], 'flooding_rp_idx_rawdata')

        # Create the output directory if it does not exist
        if not os.path.exists(output_path):

            os.makedirs(output_path, exist_ok=True)

        # Output the raster file
        with rasterio.open(
            os.path.join(output_path, output_file_name),
            'w',
            driver='GTiff',
            height=raster_array.shape[0],
            width=raster_array.shape[1],
            count=1,
            dtype=raster_array.dtype,
            crs=file_crs,
            transform=file_transform,
            compress='lzw',
        ) as f:

            # Write the raster file
            f.write(raster_array, 1)

            # Set the description of the raster file
            f.update_tags(
                1,
                NAME=variable_name,
            )
        
        return
    
    def _calc_flooding_stats(self, raster_data_multi_rp: np.ndarray) -> np.ndarray:

        """
        This method calculates the statistics of the flooding.
        ================================================================================

        Arguments:

            raster_data_multi_rp (np.ndarray): The array of the raster data with multiple return periods. The return period is the first dimension.

        Returns:

            flooding_dummy (np.ndarray): The array of the flooding dummy. The value is 1 if the flooding depth is larger than the threshold, and 0 otherwise.

            flooding_rp_idx (np.ndarray): The index of the return period that the flooding occurs. The value is 255 if the flooding depth is smaller than the threshold.
        """

        # Calculate the flooding dummy array
        flooding_dummy = np.where(raster_data_multi_rp > self.flooding_depth_threshold * 10, 1, 0).astype(np.uint8)

        # Calculate the index of the return period that the flooding occurs
        flooding_rp_idx = np.argmax(flooding_dummy, axis=0).astype(np.uint8)
        flooding_rp_idx = np.where(np.max(flooding_dummy, axis=0) == 0, 255, flooding_rp_idx)

        #return flooding_dummy, flooding_rp_idx
        return flooding_rp_idx