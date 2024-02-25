import numpy as np
import os
import itertools as it
import numba as nb
import rasterio
from . import analysis_process
from .. import environment as env

env.init()

class MergeVariant(analysis_process.AnalysisProcess):
    
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

        self.merge_variant()

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
    
    def merge_variant(self) -> None:

        """
        This method processes the merging of the variant files.
        For the files of coastal flooding (file name starts with "flooding_rp_idx.inuncoast"), the variant of subsidence condition and sea level rise scenario are merged.
        For the files of riverine flooding (file name starts with "flooding_rp_idx.inunriver"), the variant of GCM are merged.
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
        
            # Set the output file
            output_file_flooding_rp_idx = f'flooding_rp_idx.inuncoast.{climate_scenario}.{year}.tif'

            # Check if the output file exists
            if (os.path.isfile(os.path.join(env.CONFIG['path']['output'], 'flooding_rp_idx_merged', output_file_flooding_rp_idx))):
                
                continue

            # Create a array to store the raster data
            raster_data_variant = np.full((len(self.list_subsidence), len(self.list_projection), self.ref_shape[0], self.ref_shape[1]), np.nan, dtype=np.float16)

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
                
                # Set the file name
                source_file_name = 'flooding_rp_idx.inuncoast.{climatescenario}.{subsidence}.{year}.{projection}.tif'.format(
                    climatescenario=climate_scenario,
                    subsidence=subsidence,
                    year=year,
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
                raster_data_variant[self.list_subsidence.index(subsidence), self.list_projection.index(projection), ...] = raster_array.squeeze()

            # Merge the variant of data
            raster_data_variant = self._merge_variant(raster_data_variant, 'median')
            raster_data_variant = np.where(np.isnan(raster_data_variant), 255, raster_data_variant.astype(np.int8))

            # ================================================================================
            #
            # Output the raster file
            #
            # ================================================================================

            # Output the raster file
            self._output_raster_file(
                output_file_flooding_rp_idx,
                raster_data_variant,
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
        
            # Set the output file
            output_file_flooding_rp_idx = f'flooding_rp_idx.inunriver.{climate_scenario}.{year}.tif'

            # Check if the output file exists
            if (os.path.isfile(os.path.join(env.CONFIG['path']['output'], 'flooding_rp_idx_merged', output_file_flooding_rp_idx))):
                
                continue

            # Create a array to store the raster data
            raster_data_variant = np.full((len(self.list_model), self.ref_shape[0], self.ref_shape[1]), np.nan, dtype=np.float16)

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
                
                # Set the file name
                source_file_name = 'flooding_rp_idx.inunriver.{climatescenario}.{model}.{year}.tif'.format(
                    climatescenario=climate_scenario,
                    model=model,
                    year=year,
                )

                # Read the raster file
                raster_array, _, _ = self._read_raster_data(source_file_name)

                # Check if the file exists
                if (raster_array is None):

                    continue

                # Print message
                print(f'Reading the file {source_file_name}')

                # Update the raster data
                raster_data_variant[self.list_model.index(model), ...] = raster_array.squeeze()

            # Merge the variant of data
            raster_data_variant = self._merge_variant(raster_data_variant, 'median')
            raster_data_variant = np.where(np.isnan(raster_data_variant), 255, raster_data_variant.astype(np.int8))

            # ================================================================================
            #
            # Output the raster file
            #
            # ================================================================================

            # Output the raster file
            self._output_raster_file(
                output_file_flooding_rp_idx,
                raster_data_variant,
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

            source_file_name = os.listdir(os.path.join(env.CONFIG['path']['output'], 'flooding_rp_idx_rawdata'))[0]

        # ================================================================================
        #
        # Read the raster file
        #
        # ================================================================================
        
        source_file_ = os.path.join(env.CONFIG['path']['output'], 'flooding_rp_idx_rawdata', source_file_name)
        
        # Check if the file exists
        if not (os.path.isfile(source_file_)):

            print(f'The file {source_file_name} does not exist.')
            
            return None, None, None

        # Read the raster file
        with rasterio.open(source_file_) as f:

            # Get the properties of the raster file
            file_crs=f.crs
            file_transform=f.transform

            # Read the raster file
            raster_array = f.read()

        # Replace 255 to np.nan
        raster_array = np.where(raster_array == 255, np.nan, raster_array).astype(np.float16)

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

        output_path = os.path.join(env.CONFIG['path']['output'], 'flooding_rp_idx_merged')

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
    
    def _merge_variant(
            self,
            raster_data_variant: np.ndarray,
            merge_stats: str='median',
        ) -> np.ndarray:

        """
        This method processes the merge of the variant files.
        ================================================================================

        Arguments:

            raster_data_variant (np.ndarray): The array of the raster data.

            merge_stats (dict of str): The statistics of each variant to merge the variant files. The default value is {"subsidence": "median", "projection": "median", "gcm": "median"}. The key is the name of the variant, and the value is the statistics to merge the variant files. The statistics can be one of the following:

                - "median": The median of the variant files is calculated.

                - "mean": The mean of the variant files is calculated.

                - "quantile_XX": The quantile XX of the variant files is calculated. For example, "quantile_05" is the 5th percentile of the variant files.

        Returns:

            raster_data_variant (np.ndarray): The array of the raster data.
        """

        # ================================================================================
        #
        # Merge the variant of data
        #
        # ================================================================================

        @nb.jit(nopython=True, parallel=True)
        def _process_merge(data, merge_stats):
            
            data_merged = np.full((data.shape[-2], data.shape[-1]), np.nan, dtype=np.float32)

            for idx_lat in nb.prange(data.shape[-2]):

                for idx_lon in nb.prange(data.shape[-1]):

                    data_sel = data[..., idx_lat, idx_lon].flatten()

                    # Skip the process if the data is all NaN
                    if (np.all(np.isnan(data_sel))):

                        continue

                    # Merge the variant of data
                    if (merge_stats == 'median'):

                        data_merged[idx_lat, idx_lon] = np.nanmedian(data_sel)

                    elif (merge_stats == 'mean'):

                        data_merged[idx_lat, idx_lon] = np.nanmean(data_sel)

            return data_merged

        # Merge the variant of data
        raster_data_variant = _process_merge(raster_data_variant.astype(np.float32), merge_stats).astype(np.float16)

        # Squueze the dimension
        raster_data_variant = raster_data_variant.squeeze()

        # ================================================================================
        #
        # Check the result
        #
        # ================================================================================
                
        # Check if all dimensions are merged
        if (np.ndim(raster_data_variant) != 2):

            raise ValueError('The dimension of the raster data is not valid.')
        
        return raster_data_variant