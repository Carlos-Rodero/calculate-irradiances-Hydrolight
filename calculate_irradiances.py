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
import numpy as np
import sys
import time
import threading


class ProcessRadFile:
    """
    Open Lroot.txt file and extract values we need it.
    Create dataframe from content of file
    Obtain irradiances Ed, Eu and El from total radiance
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

    def open_file(self):
        """
        Open file and get content of file
        """
        f = os.path.join(self.path_files_raw, self.file_name)
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
        Calculate irradiances Ed, Eu and El
        Ed from theta values: 0 - 90, and phi values: 0 - 360
        Eu from theta values: 90 - 180, and phi values: 0 - 360
        El1 from theta values: 0 - 180, and phi values: 0 - 180
        El2 from theta values: 0 - 180, and phi values: 180 - 360

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
                    'calculated_El1': 0,
                    'calculated_El2': 0})

        df_final = pd.DataFrame(x, columns=(
            'lambda',
            'depth',
            'calculated_Ed',
            'calculated_Eu',
            'calculated_El1',
            'calculated_El2',
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
        total_Ed = 0
        total_Eu = 0
        total_El1 = 0
        total_El2 = 0
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
                total_Ed = 0
                total_Eu = 0
                total_El1 = 0
                total_El2 = 0

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

            # Calculate El1
            if ((theta >= 0) and (theta <= 180)) \
               and ((phi >= 0) and (phi <= 180)):

                # calculate in two polar cap radiance
                if (theta == 0) and (phi == 0):
                    L1_10_1 = self.df['total_radiance'].iloc[i]*2*math.pi*(
                        1-math.sin(math.radians(5)))

                if (theta == 180) and (phi == 180):
                    L1_minus10_1 = self.df['total_radiance'].iloc[i] * \
                        2*math.pi*(1-math.sin(math.radians(5)))

                # calculate over the non-polar quad. dmu from 1 to
                # minus 9, dphi from 1 to 12.

                # calculate dmu and dphi and add to cos_weight_rad
                if (theta == 87.5):
                    dmu = math.sin(math.radians(theta - 2.5)) - math.sin(
                        math.radians(theta + 2.5))
                    theta_r = theta

                    dphi = (15*2*math.pi)/360

                    sin_weight_rad_El1 += self.df['total_radiance'].iloc[i] * \
                        math.sin(math.radians(theta_r)) * dmu * dphi

                elif (theta == 92.5):
                    dmu = math.sin(math.radians(theta - 2.5)) - math.sin(
                        math.radians(theta + 2.5))
                    theta_r = theta

                    dphi = (15*2*math.pi)/360

                    sin_weight_rad_El1 += self.df['total_radiance'].iloc[i] * \
                        math.sin(math.radians(theta_r)) * dmu * dphi

                elif ((theta > 0) and (theta < 180)):
                    dmu = math.sin(math.radians(theta - 5)) - math.sin(
                        math.radians(theta + 5))
                    theta_r = theta

                    dphi = (15*2*math.pi)/360

                    sin_weight_rad_El1 += self.df['total_radiance'].iloc[i] * \
                        math.sin(math.radians(theta_r)) * dmu * dphi

                else:
                    pass

            try:
                total_El1 = sin_weight_rad_El1 + L1_10_1 + L1_minus10_1
                df_final.loc[
                    (df_final['depth'] == d) & (
                        df_final['lambda'] == lmbd),
                    'calculated_El1'] = total_El1
            except UnboundLocalError:
                pass

            # Calculate El2
            if ((theta >= 0) and (theta <= 180)) \
               and ((phi >= 180) and (phi <= 360)):

                # calculate in two polar cap radiance
                if (theta == 0) and (phi == 180):
                    L2_10_1 = self.df['total_radiance'].iloc[i]*2*math.pi*(
                        1-math.sin(math.radians(5)))

                if (theta == 180) and (phi == 180):
                    L2_minus10_1 = self.df['total_radiance'].iloc[i] * \
                        2*math.pi*(1-math.sin(math.radians(5)))

                # calculate over the non-polar quad. dmu from 1 to
                # minus 9, dphi from 1 to 12.

                # calculate dmu and dphi and add to cos_weight_rad
                if (theta == 87.5):
                    dmu = math.sin(math.radians(theta - 2.5)) - math.sin(
                        math.radians(theta + 2.5))
                    theta_r = theta

                    dphi = (15*2*math.pi)/360

                    sin_weight_rad_El2 += self.df['total_radiance'].iloc[i] * \
                        math.sin(math.radians(theta_r)) * dmu * dphi

                elif (theta == 92.5):
                    dmu = math.sin(math.radians(theta - 2.5)) - math.sin(
                        math.radians(theta + 2.5))
                    theta_r = theta

                    dphi = (15*2*math.pi)/360

                    sin_weight_rad_El2 += self.df['total_radiance'].iloc[i] * \
                        math.sin(math.radians(theta_r)) * dmu * dphi

                elif ((theta > 0) and (theta < 180)):
                    dmu = math.sin(math.radians(theta - 5)) - math.sin(
                        math.radians(theta + 5))
                    theta_r = theta

                    dphi = (15*2*math.pi)/360

                    sin_weight_rad_El2 += self.df['total_radiance'].iloc[i] * \
                        math.sin(math.radians(theta_r)) * dmu * dphi

                else:
                    pass

            try:
                total_El2 = sin_weight_rad_El2 + L2_10_1 + L2_minus10_1
                df_final.loc[
                    (df_final['depth'] == d) & (
                        df_final['lambda'] == lmbd),
                    'calculated_El2'] = total_El2
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


if __name__ == "__main__":

    prf = ProcessRadFile()
    prf.open_file()
    prf.create_dataframe_from_Lroot()
    prf.calculate_irradiances()
