# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 10:11:00 2020

@author: Carlos Rodero García

Process file Lroot from Hydrolight to obtain the irradiance (Ed, Eu or El).

"""
import os
import re
import pandas as pd
import io
import math
import sys
import time
import threading
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import numpy as np


class ProcessRadFile:
    """
    Open Lroot.txt file and extract values we need it.
    Create dataframe from content of file
    Obtain irradiances Ed, Eu and different El from total radiance
    """

    def __init__(self):

        # class variables
        self.start_string_Lroot = (r"L_dif is in-water diffuse radiance; "
                                   r"theta = 0 to 180 deg; in water only")
        self.stop_string_Lroot = r""
        self.file_name = "Lroot.txt"
        self.path_files_raw = "files/raw"
        self.path_files_csv = "files/csv"
        self.content = None
        self.df = pd.DataFrame()
        self.situation = 1

    def open_file(self, file_name=None, path_file=None):
        """
        Open file and get content of file

        Parameters
        ----------
            file_name: str
                Name of the file (Default=None)
            path_file: str
                Path of the file (Default=None)
        """
        if file_name is None:
            f = os.path.join(self.path_files_raw, self.file_name)
        else:
            f = os.path.join(path_file, file_name)
        try:
            with open(f, 'r') as file:
                self.content = file.read()

        except FileNotFoundError:
            print(f"File {self.file_name} not found")
            exit()

    def _create_dataframe_from_Lroot(self):
        """
        Create dataframe from content file
        """

        # Select the useful data from content. \Z to select EOF in RE
        patron = re.compile(r'{}(?P<table>[\s\S]*?){}\Z'.format(
                            self.start_string_Lroot,
                            self.stop_string_Lroot))
        selected_info = r""

        for m in re.finditer(patron, self.content):
            selected_info = m.group('table')

        # Replace column name
        selected_info = re.sub('total radiance L_sky or L_dir',
                               'total radiance, L_sky or L_dir',
                               selected_info)

        col_names = ["depth", "theta", "phi", "lambda", "total_radiance",
                     "L_sky-or-L_dir", "L_w-or-L_dif", "L_sr"]

        self.df = pd.read_fwf(io.StringIO(selected_info), header=1,
                              skipinitialspace=True)

        self.df.columns = col_names

        # Delete unused rows
        self.df.drop(index=0, inplace=True)

        # Save as csv
        fname = f"{self.file_name.split('.')[0]}_data.csv"
        f = os.path.join(self.path_files_csv, fname)
        self.df.to_csv(f)

    def create_dataframe_from_Lroot(self):
        """
        Create new thread to process create dataframe
        """
        the_process = threading.Thread(
            target=self._create_dataframe_from_Lroot)
        the_process.start()
        while the_process.is_alive():
            self.animated_loading(process_name="Create Dataframe")

    def _calculate_irradiances(self):
        """
        Calculate irradiances Ed, Eu and different El
        - Ed from theta values: 0 - 90, and phi values: 0 - 360
        - Eu from theta values: 90 - 180, and phi values: 0 - 360
        - El1 without polar cap from theta values: 0 - 180,
        and phi values: 0 - 180
        - El2 without polar capfrom theta values: 0 - 180,
        and phi values: 180 - 360
        - El1 with polar cap from theta values: 0 - 180,
        and phi values: 0 - 180
        - El2 with polar cap from theta values: 0 - 180,
        and phi values: 180 - 360
        - Ehc (horizontal crown) from theta values: 87.5 - 92.5,
        and phi values: 0 - 360
        - Ehc45 (horizontal crown at 45º) from theta values: 40 - 50,
        and phi values: 0 - 360

        Return
        ------
            df_final: pandas dataframe object
                dataframe with Ed values
        """
        # Creation of pandas dataframe with useful data
        x = []
        for i in self.df['lambda'].unique():
            for z in self.df['depth'].unique():
                x.append({
                    'lambda': i,
                    'depth': z,
                    'calculated_Ed': 0,
                    'calculated_Eu': 0,
                    'calculated_El1_no_polar_cap': 0,
                    'calculated_El2_no_polar_cap': 0,
                    'calculated_El1_polar_cap': 0,
                    'calculated_El2_polar_cap': 0,
                    'calculated_Ehc': 0,
                    'calculated_Ehc_45': 0})

        df_final = pd.DataFrame(x, columns=(
            'lambda',
            'depth',
            'calculated_Ed',
            'calculated_Eu',
            'calculated_El1_no_polar_cap',
            'calculated_El2_no_polar_cap',
            'calculated_El1_polar_cap',
            'calculated_El2_polar_cap',
            'calculated_Ehc',
            'calculated_Ehc_45'
            ))

        df_final['lambda'] = df_final['lambda'].astype(float).fillna(0.0)
        df_final['depth'] = df_final['depth'].astype(float).fillna(0.0)

        self.df = self.df.apply(pd.to_numeric, args=('coerce',))

        # declare variables
        # declare first depth in dataframe
        d = -1
        lmbd = self.df['lambda'].iloc[0]
        cos_weight_rad_Ed = 0
        cos_weight_rad_Eu = 0
        sin_weight_rad_El1 = 0
        sin_weight_rad_El2 = 0
        sin_weight_rad_El1_polar_cap = 0
        sin_weight_rad_El2_polar_cap = 0
        sin_weight_rad_Ehc = 0
        sin_weight_rad_Ehc_45 = 0
        total_Ed = 0
        total_Eu = 0
        total_El1 = 0
        total_El2 = 0
        total_El1_polar_cap = 0
        total_El2_polar_cap = 0
        total_Ehc = 0
        total_Ehc_45 = 0
        dmu = 0
        dphi = 0
        theta_r = 0

        # Calculate irradiance for every depth and every lambda
        print(f" - Start in lambda = {lmbd}")

        for i in range(0, len(self.df)):
            theta = self.df['theta'].iloc[i]
            phi = self.df['phi'].iloc[i]
            depth = self.df['depth'].iloc[i]

            # clear variables in each new depth
            if d != depth:
                d = self.df['depth'].iloc[i]
                cos_weight_rad_Ed = 0
                cos_weight_rad_Eu = 0
                sin_weight_rad_El1 = 0
                sin_weight_rad_El2 = 0
                sin_weight_rad_El1_polar_cap = 0
                sin_weight_rad_El2_polar_cap = 0
                sin_weight_rad_Ehc = 0
                sin_weight_rad_Ehc_45 = 0
                total_Ed = 0
                total_Eu = 0
                total_El1 = 0
                total_El2 = 0
                total_El1_polar_cap = 0
                total_El2_polar_cap = 0
                total_Ehc = 0
                total_Ehc_45 = 0

            if lmbd != self.df['lambda'].iloc[i]:
                lmbd = self.df['lambda'].iloc[i]
                print(f" - Calculate in lambda = {lmbd}")

            # Calculate Ed
            if ((theta >= 0) and (theta <= 90)) \
               and ((phi >= 0) and (phi <= 360)):

                # calculate in polar cap radiance
                if (theta == 0) and (phi == 0):
                    L_10_1 = self.df['total_radiance'].iloc[i]*2*math.pi*(
                        1-math.cos(math.radians(5)))

                # calculate over the non-polar quad. dmu from 1 to 9,
                # dphi from 1 to 24.

                # calculate dmu and dphi and add to cos_weight_rad
                if (theta == 87.5):
                    dmu = math.cos(math.radians(theta - 2.5)) - math.cos(
                        math.radians(theta + 2.5))
                    theta_r = theta

                    dphi = (15*2*math.pi)/360

                    cos_weight_rad_Ed += self.df['total_radiance'].iloc[i] * \
                        math.cos(math.radians(theta_r)) * dmu * dphi

                elif ((theta > 0) and (theta < 87.5)):
                    dmu = math.cos(math.radians(theta - 5)) - math.cos(
                        math.radians(theta + 5))
                    theta_r = theta

                    dphi = (15*2*math.pi)/360

                    cos_weight_rad_Ed += self.df['total_radiance'].iloc[i] * \
                        math.cos(math.radians(theta_r)) * dmu * dphi

                else:
                    pass

            try:
                total_Ed = cos_weight_rad_Ed + L_10_1
                df_final.loc[
                    (df_final['depth'] == d) & (
                        df_final['lambda'] == lmbd),
                    'calculated_Ed'] = total_Ed
            except UnboundLocalError:
                pass

            # Calculate Eu
            if ((theta >= 90) and (theta <= 180)) \
               and ((phi >= 0) and (phi <= 360)):

                # calculate in polar cap radiance
                if (theta == 180) and (phi == 180):
                    L_minus10_1 = self.df['total_radiance'].iloc[i]*2*math.pi*(
                        1-math.cos(math.radians(5)))

                # calculate over the non-polar quad. dmu from minus 1 to
                # minus 9, dphi from 1 to 24.

                # calculate dmu and dphi and add to cos_weight_rad
                if (theta == 92.5):
                    dmu = math.cos(math.radians(theta - 2.5)) - math.cos(
                        math.radians(theta + 2.5))
                    theta_r = theta

                    dphi = (15*2*math.pi)/360

                    cos_weight_rad_Eu += self.df['total_radiance'].iloc[i] * \
                        math.cos(math.radians(theta_r)) * dmu * dphi

                elif ((theta > 90) and (theta < 180)):
                    dmu = math.cos(math.radians(theta - 5)) - math.cos(
                        math.radians(theta + 5))
                    theta_r = theta

                    dphi = (15*2*math.pi)/360

                    cos_weight_rad_Eu += self.df['total_radiance'].iloc[i] * \
                        math.cos(math.radians(theta_r)) * dmu * dphi

                else:
                    pass

            try:
                total_Eu = cos_weight_rad_Eu + L_minus10_1
                df_final.loc[
                    (df_final['depth'] == d) & (
                        df_final['lambda'] == lmbd),
                    'calculated_Eu'] = -1 * total_Eu
            except UnboundLocalError:
                pass

            # Situation 1. El1 hemisphere: -90 < phi < 90
            if self.situation == 1:

                # Calculate El1_no_polar_cap and El1_polar_cap

                if ((theta >= 0) and (theta <= 180)) \
                   and ((phi > 0) and (phi <= 90)) or \
                       ((phi >= 270) and (phi <= 360)):
                    # In phi = 0 or phi = 360 we have the same radiance
                    # measurement, for this reason we only take one

                    # calculate in two polar cap radiance.
                    # Divide by 2 to get half polar cap radiance
                    if (theta == 0) and (phi == 15):
                        L1_10_1 = self.df['total_radiance'].iloc[i]*2*math.pi*(
                            1-math.cos(math.radians(5)))/2

                    elif (theta == 180) and (phi == 15):
                        L1_minus10_1 = (
                            self.df['total_radiance'].iloc[i]*2*math.pi*(
                                1-math.cos(math.radians(5)))/2)

                    # calculate over the non-polar quad. dmu from 1 to
                    # minus 9, dphi from 1 to 12.

                    # calculate dmu and dphi and add to sin_weight_rad_El1
                    elif ((theta == 87.5) or (theta == 92.5)):
                        dmu = math.cos(math.radians(theta - 2.5)) - math.cos(
                            math.radians(theta + 2.5))
                        theta_r = theta
                        phi_s = phi

                        dphi = (15*2*math.pi)/360

                        sin_weight_rad_El1 += self.df[
                            'total_radiance'].iloc[i] * \
                            abs(math.sin(math.radians(theta_r)) *
                                math.cos(math.radians(phi_s))) * dmu * dphi

                    elif ((theta > 0) and (theta < 180)):
                        dmu = math.cos(math.radians(theta - 5)) - math.cos(
                            math.radians(theta + 5))
                        theta_r = theta
                        phi_s = phi

                        dphi = (15*2*math.pi)/360

                        sin_weight_rad_El1 += self.df[
                            'total_radiance'].iloc[i] * \
                            abs(math.sin(math.radians(theta_r)) *
                                math.cos(math.radians(phi_s))) * dmu * dphi

                    else:
                        pass

                try:
                    # with polar cap
                    sin_weight_rad_El1_polar_cap = (
                        sin_weight_rad_El1 + L1_10_1 + L1_minus10_1)

                    total_El1_polar_cap = sin_weight_rad_El1_polar_cap

                    df_final.loc[
                        (df_final['depth'] == d) & (
                            df_final['lambda'] == lmbd),
                        'calculated_El1_polar_cap'] = total_El1_polar_cap

                    # without polar cap
                    total_El1 = sin_weight_rad_El1

                    df_final.loc[
                        (df_final['depth'] == d) & (
                            df_final['lambda'] == lmbd),
                        'calculated_El1_no_polar_cap'] = total_El1
                except UnboundLocalError:
                    pass

                # Calculate El2_no_polar_cap and El2_polar_cap

                if ((theta >= 0) and (theta <= 180)) \
                   and ((phi >= 90) and (phi <= 270)):

                    # calculate in two polar cap radiance
                    # Divide by 2 to get half polar cap radiance
                    if (theta == 0) and (phi == 180):
                        L2_10_1 = self.df['total_radiance'].iloc[i]*2*math.pi*(
                            1-math.cos(math.radians(5)))/2

                    elif (theta == 180) and (phi == 180):
                        L2_minus10_1 = (
                            self.df['total_radiance'].iloc[i]*2*math.pi*(
                                1-math.cos(math.radians(5)))/2)

                    # calculate over the non-polar quad. dmu from 1 to
                    # minus 9, dphi from 1 to 12.

                    # calculate dmu and dphi and add to sin_weight_rad_El2
                    elif ((theta == 87.5) or (theta == 92.5)):
                        dmu = math.cos(math.radians(theta - 2.5)) - math.cos(
                            math.radians(theta + 2.5))
                        theta_r = theta
                        phi_s = phi

                        dphi = (15*2*math.pi)/360

                        sin_weight_rad_El2 += self.df[
                            'total_radiance'].iloc[i] * \
                            abs(math.sin(math.radians(theta_r)) *
                                math.cos(math.radians(phi_s))) * dmu * dphi

                    elif ((theta > 0) and (theta < 180)):
                        dmu = math.cos(math.radians(theta - 5)) - math.cos(
                            math.radians(theta + 5))
                        theta_r = theta
                        phi_s = phi

                        dphi = (15*2*math.pi)/360

                        sin_weight_rad_El2 += self.df[
                            'total_radiance'].iloc[i] * \
                            abs(math.sin(math.radians(theta_r)) *
                                math.cos(math.radians(phi_s))) * dmu * dphi

                    else:
                        pass

                try:
                    # with polar cap
                    sin_weight_rad_El2_polar_cap = (
                        sin_weight_rad_El2 + L2_10_1 + L2_minus10_1)

                    total_El2_polar_cap = sin_weight_rad_El2_polar_cap

                    df_final.loc[
                        (df_final['depth'] == d) & (
                            df_final['lambda'] == lmbd),
                        'calculated_El2_polar_cap'] = total_El2_polar_cap

                    # without polar cap
                    total_El2 = sin_weight_rad_El2

                    df_final.loc[
                        (df_final['depth'] == d) & (
                            df_final['lambda'] == lmbd),
                        'calculated_El2_no_polar_cap'] = total_El2
                except UnboundLocalError:
                    pass

            # Situation 2. El1 hemisphere: 0 < phi < 180
            elif self.situation == 2:
                # Calculate El1_no_polar_cap and El1_polar_cap

                if ((theta >= 0) and (theta <= 180)) \
                   and ((phi >= 0) and (phi <= 180)):

                    # calculate in two polar cap radiance
                    # Divide by 2 to get half polar cap radiance
                    if (theta == 0) and (phi == 0):
                        L1_10_1 = (
                            self.df['total_radiance'].iloc[i]*2*math.pi*(
                                1-math.sin(math.radians(5))))/2

                    if (theta == 180) and (phi == 180):
                        L1_minus10_1 = (
                            self.df['total_radiance'].iloc[i] *
                            2*math.pi*(1-math.sin(math.radians(5))))/2

                    # calculate over the non-polar quad. dmu from 1 to
                    # minus 9, dphi from 1 to 12.

                    # calculate dmu and dphi and add to sin_weight_rad_El1
                    if ((theta == 87.5) or (theta == 92.5)):
                        dmu = math.cos(math.radians(theta - 2.5)) - math.cos(
                            math.radians(theta + 2.5))
                        theta_r = theta
                        phi_s = phi

                        dphi = (15*2*math.pi)/360

                        sin_weight_rad_El1 += self.df[
                            'total_radiance'].iloc[i] * \
                            abs(math.sin(math.radians(theta_r)) *
                                math.sin(math.radians(phi_s))) * dmu * dphi

                    elif ((theta > 0) and (theta < 180)):
                        dmu = math.cos(math.radians(theta - 5)) - math.cos(
                            math.radians(theta + 5))
                        theta_r = theta
                        phi_s = phi

                        dphi = (15*2*math.pi)/360

                        sin_weight_rad_El1 += self.df[
                            'total_radiance'].iloc[i] * \
                            abs(math.sin(math.radians(theta_r)) *
                                math.sin(math.radians(phi_s))) * dmu * dphi

                    else:
                        pass

                try:
                    # with polar cap
                    sin_weight_rad_El1_polar_cap = (
                        sin_weight_rad_El1 + L1_10_1 + L1_minus10_1)

                    total_El1_polar_cap = sin_weight_rad_El1_polar_cap

                    df_final.loc[
                        (df_final['depth'] == d) & (
                            df_final['lambda'] == lmbd),
                        'calculated_El1_polar_cap'] = total_El1_polar_cap

                    # without polar cap
                    total_El1 = sin_weight_rad_El1

                    df_final.loc[
                        (df_final['depth'] == d) & (
                            df_final['lambda'] == lmbd),
                        'calculated_El1_no_polar_cap'] = total_El1
                except UnboundLocalError:
                    pass

                # Calculate El2_no_polar_cap and El2_polar_cap

                if ((theta >= 0) and (theta <= 180)) \
                   and ((phi >= 180) and (phi <= 360)):

                    # calculate in two polar cap radiance
                    # Divide by 2 to get half polar cap radiance
                    if (theta == 0) and (phi == 180):
                        L2_10_1 = (
                            self.df['total_radiance'].iloc[i]*2*math.pi*(
                                1-math.sin(math.radians(5))))/2

                    if (theta == 180) and (phi == 180):
                        L2_minus10_1 = (
                            self.df['total_radiance'].iloc[i] *
                            2*math.pi*(1-math.sin(math.radians(5))))/2

                    # calculate over the non-polar quad. dmu from 1 to
                    # minus 9, dphi from 1 to 12.

                    # calculate dmu and dphi and add to sin_weight_rad_El2
                    if ((theta == 87.5) or (theta == 92.5)):
                        dmu = math.cos(math.radians(theta - 2.5)) - math.cos(
                            math.radians(theta + 2.5))
                        theta_r = theta
                        phi_s = phi

                        dphi = (15*2*math.pi)/360

                        sin_weight_rad_El2 += self.df[
                            'total_radiance'].iloc[i] * \
                            abs(math.sin(math.radians(theta_r)) *
                                math.sin(math.radians(phi_s))) * dmu * dphi

                    elif ((theta > 0) and (theta < 180)):
                        dmu = math.cos(math.radians(theta - 5)) - math.cos(
                            math.radians(theta + 5))
                        theta_r = theta
                        phi_s = phi

                        dphi = (15*2*math.pi)/360

                        sin_weight_rad_El2 += self.df[
                            'total_radiance'].iloc[i] * \
                            abs(math.sin(math.radians(theta_r)) *
                                math.sin(math.radians(phi_s))) * dmu * dphi

                    else:
                        pass

                try:
                    # with polar cap
                    sin_weight_rad_El2_polar_cap = (
                        sin_weight_rad_El2 + L2_10_1 + L2_minus10_1)

                    total_El2_polar_cap = sin_weight_rad_El2_polar_cap

                    df_final.loc[
                        (df_final['depth'] == d) & (
                            df_final['lambda'] == lmbd),
                        'calculated_El2_polar_cap'] = total_El2_polar_cap

                    # without polar cap
                    total_El2 = sin_weight_rad_El2  # + L2_10_1 + L2_minus10_1

                    df_final.loc[
                        (df_final['depth'] == d) & (
                            df_final['lambda'] == lmbd),
                        'calculated_El2_no_polar_cap'] = total_El2
                except UnboundLocalError:
                    pass

            # Calculate Ehc horizontal crown
            if ((theta >= 87.5) and (theta <= 92.5)) \
               and ((phi >= 0) and (phi <= 360)):

                # calculate over the non-polar quad. dmu from 1 to -1,
                # dphi from 1 to 24.

                # calculate dmu and dphi and add to sin_weight_rad_Ehc
                if ((theta == 87.5) or (theta == 92.5)):
                    dmu = math.cos(math.radians(theta - 2.5)) - math.cos(
                        math.radians(theta + 2.5))
                    theta_r = theta
                    phi_s = phi

                    dphi = (15*2*math.pi)/360

                    sin_weight_rad_Ehc += self.df[
                        'total_radiance'].iloc[i] * \
                        abs(math.sin(math.radians(theta_r)) *
                            math.cos(math.radians(phi_s))) * dmu * dphi

                else:
                    pass

            try:
                total_Ehc = sin_weight_rad_Ehc
                df_final.loc[
                    (df_final['depth'] == d) & (
                        df_final['lambda'] == lmbd),
                    'calculated_Ehc'] = total_Ehc
            except UnboundLocalError:
                pass

            # Calculate Ehc horizontal crown 45º
            if ((theta >= 40) and (theta <= 50)) \
               and ((phi >= 0) and (phi <= 360)):

                # calculate over the non-polar quad. dmu from 4 to 5,
                # dphi from 1 to 24.

                if ((theta >= 40) and (theta <= 50)):
                    dmu = math.cos(math.radians(theta - 5)) - math.cos(
                        math.radians(theta + 5))
                    theta_r = theta
                    phi_s = phi

                    dphi = (15*2*math.pi)/360

                    sin_weight_rad_Ehc_45 += self.df[
                        'total_radiance'].iloc[i] * \
                        abs(math.sin(math.radians(theta_r)) *
                            math.cos(math.radians(phi_s))) * dmu * dphi

                else:
                    pass

            try:
                total_Ehc_45 = sin_weight_rad_Ehc_45
                df_final.loc[
                    (df_final['depth'] == d) & (
                        df_final['lambda'] == lmbd),
                    'calculated_Ehc_45'] = total_Ehc_45
            except UnboundLocalError:
                pass

        # save as csv
        fname = f"{self.file_name.split('.')[0]}_calculated_irradiances.csv"
        f = os.path.join(self.path_files_csv, fname)
        df_final.to_csv(f)

        return df_final

    def calculate_irradiances(self):
        the_process = threading.Thread(
            target=self._calculate_irradiances)
        start = time.time()
        the_process.start()
        while the_process.is_alive():
            self.animated_loading(process_name="Calculate irradiances")
        end = time.time()
        print("\nComplete. ")
        print(f"Time calculating irradiances: {(end - start)/60} minutes")

    def animated_loading(self, process_name=""):
        chars = r"/—\|"
        for char in chars:
            sys.stdout.write('\r'+f'{process_name} - Loading...'+char)
            time.sleep(.1)
            sys.stdout.flush()

    def _create_dataframe_from_Lroot_data(self):
        """
        Create dataframe from content file
        """
        # Create dataframe from content of .csv file
        self.df = pd.read_csv(io.StringIO(self.content), header=0,
                              skipinitialspace=True, index_col=0)

    def plot_heatmap_radiances(self,  has_show=False):
        """
        Plot radiances as heatmaps in different depths and lambdas in Plotly

        Parameters
        ----------
            has_show: Boolean
                Flag to show the plot. By default, False.

        """
        self.df = self.df.apply(pd.to_numeric, args=('coerce',))
        df = self.df
        list_radiance = []
        index = None
        index_row = None
        index_column = None

        fig = make_subplots(
            rows=2, cols=3,
            column_widths=[0.5, 0.5, 0.5],
            row_heights=[0.5, 0.5],
            subplot_titles=("Radiances at 460",
                            "Radiances at 500",
                            "Radiances at 540",
                            "Radiances at 560",
                            "Radiances at 600",
                            "Radiances at 640"))

        # create heatmap with following axes (to represent the hydrolight
        # particle in 2D)
        # x = phi   (25 elements)
        # y = theta (20 elements)
        # z = list of lists of radiance for each theta (there are 25 elements
        # for each list)
        for z in df['depth'].unique():
            if (z == 0.25) or (z == 1.00) or (z == 5.00):
                index_row = 1
                index_column = 1
                index = 1

                for i in df['lambda'].unique():
                    if (i == 460) or (i == 500) or (i == 540) or \
                     (i == 560) or (i == 600) or (i == 640):
                        # clear data
                        list_radiance.clear()
                        # count
                        if index_row == 3:
                            index_row = 1
                        if index_column == 4:
                            index_column = 1

                        for t in df['theta'].unique():

                            radiance = list(
                                df['total_radiance'].loc[(
                                    df['lambda'] == i) & (
                                        df['depth'] == z) & (
                                            df['theta'] == t)])
                            list_radiance.append(radiance)

                        # fig = go.Figure(data=go.Heatmap(
                        data = go.Heatmap(
                            z=list_radiance, zmin=0, zmax=1,
                            x=df['phi'].unique(),
                            x0=0,
                            dx=15,
                            y=df['theta'].unique(),
                            y0=0,
                            dy=10,
                            colorbar=dict(title='Range'),
                            hovertemplate='phi: %{x}<br>theta: %{y}' +
                            '<br>radiance: %{z}<extra></extra>',
                            hoverongaps=False)

                        fig.append_trace(data, index_row, index_column)

                        fig.update_xaxes(
                            title=go.layout.xaxis.Title(
                                text='phi'))

                        fig.update_yaxes(
                            autorange="reversed",
                            title=go.layout.yaxis.Title(
                                text='theta'))

                        fig.update_layout(
                            title=f'Radiances at depth: {z}')

                        index_row += 1
                        index_column += 1
                        index += 1

                if has_show is True:
                    fig.show(config={'showLink': True})

                if not os.path.exists("images/plotly"):
                    os.mkdir("images/plotly")
                fig.write_image(
                    f"images/plotly/heatmap_radiance_depth_{z}.svg")
                fig.write_html(
                    f"images/plotly/heatmap_radiance_depth_{z}.html")

        # to do
        # heatmap as a circular particle of Hydrolight. Example below
        """ N = 300
        R = 1
        x = np.linspace(-R, R, N)
        y = np.linspace(-R, R, N)
        X, Y = np.meshgrid(x, y)

        disk = X**2+Y**2
        I, J = np.where(disk > R)
        z = X*Y**2-X**3*np.sin(X-Y)
        # mask the outside of the disk of center (0,0)   and radius R
        z[I, J] = None

        trace = dict(type='heatmap',
                     x=x,
                     y=y,
                     # note that z has the shape of X,Y,
                     # not x,y as in your exmple!
                     z=z,
                     colorscale='YlGnBu',
                     showscale=True,
                     colorbar=dict(thickness=20,
                                   ticklen=4,
                                   tickfont=dict(size=10)))

        layout = dict(title='Polar heatmap',
                      width=450,
                      height=450,
                      showlegend=False,
                      xaxis=dict(visible=False),
                      yaxis=dict(visible=False)
                      )

        fig = go.Figure(layout=layout)
        fig.add_trace(trace)
        fig.write_html("example.html") """


if __name__ == "__main__":

    prf = ProcessRadFile()
    # prf.open_file()
    # prf.create_dataframe_from_Lroot()
    # prf.calculate_irradiances()

    # Flag to show or hide plot
    has_show = True

    # plot radiances in a heatmap
    prf.open_file(
        file_name="Lroot_data.csv",
        path_file="files/csv")
    prf._create_dataframe_from_Lroot_data()
    prf.plot_heatmap_radiances()
