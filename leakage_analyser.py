import h5py
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
from scipy import stats as stats
from tqdm import tqdm
import math
import os

class trace_analyser:
    def __init__(self):
        self.name = None
        self.file = None
        self.keys = []
        self.x = 0
        self.plain = None
        self.nb_traces = 0
        self.trace_len = 0
        self.traces = None
        self.target_byte = 0
        self.plain_len = 0
        self.file_x = None

        #Option settings
        self.draw_diagram = True
        self.generate_csv = True
        self.return_results = False

        #Predefined Sbox
        self.AES_Sbox = np.array([
            0x63, 0x7C, 0x77, 0x7B, 0xF2, 0x6B, 0x6F, 0xC5, 0x30, 0x01, 0x67, 0x2B, 0xFE, 0xD7, 0xAB, 0x76,
            0xCA, 0x82, 0xC9, 0x7D, 0xFA, 0x59, 0x47, 0xF0, 0xAD, 0xD4, 0xA2, 0xAF, 0x9C, 0xA4, 0x72, 0xC0,
            0xB7, 0xFD, 0x93, 0x26, 0x36, 0x3F, 0xF7, 0xCC, 0x34, 0xA5, 0xE5, 0xF1, 0x71, 0xD8, 0x31, 0x15,
            0x04, 0xC7, 0x23, 0xC3, 0x18, 0x96, 0x05, 0x9A, 0x07, 0x12, 0x80, 0xE2, 0xEB, 0x27, 0xB2, 0x75,
            0x09, 0x83, 0x2C, 0x1A, 0x1B, 0x6E, 0x5A, 0xA0, 0x52, 0x3B, 0xD6, 0xB3, 0x29, 0xE3, 0x2F, 0x84,
            0x53, 0xD1, 0x00, 0xED, 0x20, 0xFC, 0xB1, 0x5B, 0x6A, 0xCB, 0xBE, 0x39, 0x4A, 0x4C, 0x58, 0xCF,
            0xD0, 0xEF, 0xAA, 0xFB, 0x43, 0x4D, 0x33, 0x85, 0x45, 0xF9, 0x02, 0x7F, 0x50, 0x3C, 0x9F, 0xA8,
            0x51, 0xA3, 0x40, 0x8F, 0x92, 0x9D, 0x38, 0xF5, 0xBC, 0xB6, 0xDA, 0x21, 0x10, 0xFF, 0xF3, 0xD2,
            0xCD, 0x0C, 0x13, 0xEC, 0x5F, 0x97, 0x44, 0x17, 0xC4, 0xA7, 0x7E, 0x3D, 0x64, 0x5D, 0x19, 0x73,
            0x60, 0x81, 0x4F, 0xDC, 0x22, 0x2A, 0x90, 0x88, 0x46, 0xEE, 0xB8, 0x14, 0xDE, 0x5E, 0x0B, 0xDB,
            0xE0, 0x32, 0x3A, 0x0A, 0x49, 0x06, 0x24, 0x5C, 0xC2, 0xD3, 0xAC, 0x62, 0x91, 0x95, 0xE4, 0x79,
            0xE7, 0xC8, 0x37, 0x6D, 0x8D, 0xD5, 0x4E, 0xA9, 0x6C, 0x56, 0xF4, 0xEA, 0x65, 0x7A, 0xAE, 0x08,
            0xBA, 0x78, 0x25, 0x2E, 0x1C, 0xA6, 0xB4, 0xC6, 0xE8, 0xDD, 0x74, 0x1F, 0x4B, 0xBD, 0x8B, 0x8A,
            0x70, 0x3E, 0xB5, 0x66, 0x48, 0x03, 0xF6, 0x0E, 0x61, 0x35, 0x57, 0xB9, 0x86, 0xC1, 0x1D, 0x9E,
            0xE1, 0xF8, 0x98, 0x11, 0x69, 0xD9, 0x8E, 0x94, 0x9B, 0x1E, 0x87, 0xE9, 0xCE, 0x55, 0x28, 0xDF,
            0x8C, 0xA1, 0x89, 0x0D, 0xBF, 0xE6, 0x42, 0x68, 0x41, 0x99, 0x2D, 0x0F, 0xB0, 0x54, 0xBB, 0x16
        ])
    #This method extracts all Sbox values from a given trace and byte
    def get_sbox_values(self):
        input_data = self.file_x["metadata"][:self.nb_traces]
        plains = input_data["plaintext"][:, :]
        keys = input_data["key"][:, :]
        plains = plains[range(self.nb_traces), self.target_byte]

        key_byte = self.target_byte % len(keys[0])
        keys = keys[range(self.nb_traces), key_byte]
        keys = np.array(keys, dtype=int)
        plains = np.array(plains, dtype=int)
        sbox_values = self.AES_Sbox[plains ^ keys]

        return sbox_values

    #This method creates the contigency table for the given trace point
    def make_contingency_table(self, trace_point):
        sbox_values = self.get_sbox_values()
        traces = self.file_x["traces"][:, :]

        # Extract the trace data for the specified trace point
        trace = traces[range(self.nb_traces), trace_point - 1]

        # Determine the minimum and maximum data values in the trace
        data_min = np.min(trace)
        data_max = np.max(trace)

        # Determine the number of data value groups to use in the contingency table
        if (data_min < 0):
            group_amount = (int)(data_max - data_min)
        else:
            group_amount = (int)(data_max + data_min)

        interval_width = (data_max - data_min) / group_amount

        # Initialize the contingency table
        contingency_table = np.zeros((256, group_amount), dtype=int)

        # Fill in the contingency table with counts of each S-box value in each data value group
        for counter, sbox in enumerate(sbox_values):
            ind = np.where(self.AES_Sbox == sbox)[0][0]
            interval_index = int((trace[counter] - data_min) / interval_width)
            interval_index = max(0, min(interval_index, group_amount - 1))
            contingency_table[ind][interval_index] = contingency_table[ind][interval_index] + 1

        group_sums = np.sum(contingency_table, axis=0)
        deleted = 0
        for counter, sum in enumerate(group_sums):

            if sum == 0:
                contingency_table = np.delete(contingency_table, counter - deleted, 1)
                deleted = deleted + 1

        return contingency_table

    def start_t_test(self):

        sbox_values = self.get_sbox_values()
        # Get the traces from the input file
        traces = self.file_x["traces"][:, :]

        # Initialize a list of RunningTtest objects for each S-box value
        t_tests = [RunningTtest() for b in range(256)]

        # Update the RunningTtest objects with the trace data
        for b in tqdm(range(256)):
            for count, trace in enumerate(traces):
                # Update the RunningTtest object corresponding to the current S-box value with the current trace data
                if (sbox_values[count] == b):
                    t_tests[b].update(trace, 0)
                else:
                    t_tests[b].update(trace, 1)

        # Compute the t-test results for each S-box value
        Ttest_results = np.array([t_tests[b]() for b in range(256)])

        # Determine if leakage occurred for each S-box value based on the t-test results
        leakage = []
        for count, tt in enumerate(Ttest_results):
            leakage.append(0)
            for t in tt:
                if t > 4.5 or t < -4.5:
                    leakage[count] = leakage[count] + 1
        return Ttest_results, leakage

    #This is the core of the chi2 test
    def start_chi2_test(self):
        # Initialize lists to store the test results
        chi2 = []
        p = []
        dof = []
        expected_values = []
        leakage = []

        # Set the significance level for the test
        alfa = 10 ** -5

        # Create a progress bar to track the test progress
        t = tqdm(total=self.trace_len)

        # Perform the chi-squared test on each trace point
        for i in range(self.trace_len):

            t.update()
            # Create a contingency table for the current trace point
            con_table = self.make_contingency_table(i)

            # Compute the chi-squared test statistic, p-value, degrees of freedom, and expected values
            chi2_i, p_i, dof_i, expected_values_i = scipy.stats.chi2_contingency(con_table)

            # Append the test results to the corresponding lists
            chi2.append(chi2_i)
            p.append(p_i)
            dof.append(dof_i)
            expected_values.append(expected_values_i)

            # Determines if leakage occurred based on the p-value and the significance level
            leakage_i = False
            if (p_i <= alfa):
                leakage_i = True
            leakage.append(leakage_i)

        return chi2, p, dof, expected_values, leakage

    #Calculates negative 10 based logarithm
    def neg_log(self, value):
        return -math.log(value, 10)

    # This method makes the diagram for the t-test
    def t_diagrams(self, result,leakage,target_text):
        for res in result:
            plt.plot(res)
        plt.xlabel("measurement times")
        plt.ylabel("t-value")
        plot_name = self.name.replace(".h5", "")
        plot_name = plot_name + "_" + self.keys[self.x] + "_" + target_text + "_Welch's t-test by measurement times"
        plt.title(plot_name)
        save_name = plot_name + ".png"
        plt.savefig(save_name)
        plt.clf()

        plt.plot(result)
        plt.xlabel("sbox values")
        plt.ylabel("t-values")
        plot_name = self.name.replace(".h5", "")
        plot_name = plot_name + "_" + self.keys[self.x] + "_" + target_text + "_Welch's t-test"
        plt.title(plot_name)
        save_name = plot_name + ".png"
        plt.savefig(save_name)
        plt.clf()

    #This method makes the diagrams for the chi2 test
    def chi_diagrams(self,p,target_text,neg_log_v):
        plt.plot(p)
        plot_name = self.name.replace(".h5", "")
        plot_name = plot_name + "_" + self.keys[self.x] + "_" + target_text + "_Pearsons-Chi2-test"
        plt.title(plot_name)
        plt.xlabel("measurement times ")
        plt.ylabel("p-values")
        save_name = plot_name + ".png"
        plt.savefig(save_name)
        plt.clf()

        plt.title(plot_name)
        lneg = neg_log_v(p)
        plt.plot(lneg)
        plt.xlabel("measurement times")
        plt.ylabel("-log10(p)")
        save_name = plot_name + "_log.png"
        plt.savefig(save_name)
        plt.clf()






    #This method takes an analysis option opt and a the target byte b, but if b isnot given
    # it will ask the user for the byte specification then executes the tests
    def option_choosed(self, opt, b=-1):
        plt.clf()
        neg_log_v = np.vectorize(self.neg_log)
        target_byte = b
        if (opt != 3 and target_byte < 0):
            while (True):
                print(f"Choose target byte between 0 and {self.plain_len}")
                try:
                    target_byte = int(input())
                    if (self.x > self.plain_len - 1 or self.x < 0):
                        continue
                    break
                except:
                    continue
        self.target_byte = target_byte
        f = plt.figure()
        f.set_figwidth(20)
        f.set_figheight(10)

        if opt == 0:
            self.file_x = self.file[f"{self.keys[self.x]}"]
            print("Starting Welch's t-test")
            sv_name = self.name.replace(".h5", "")
            result, leakage = self.start_t_test()
            target_text = str(target_byte)
            if self.generate_csv:
                np.savetxt(f"{sv_name}_{ self.keys[self.x]}_{target_text}_t-test_t.csv", result)

            print(f'Leakage was detected in {sum(leakage)} points')
            self.t_diagrams(result, leakage, target_text)
            if self.return_results:
                return result, leakage, None, None, None, None, None

        elif opt == 1:
            self.file_x = self.file[f"{self.keys[self.x]}"]
            chi2, p, dof, expected_values, leakage2 = self.start_chi2_test()
            target_text = str(target_byte)
            sv_name = self.name.replace(".h5", "")

            if self.generate_csv:
                np.savetxt(f"{sv_name}_{target_text}_Chi2.csv", chi2)
                np.savetxt(f"{sv_name}_{target_text}_Chi2_p.csv", p)
                np.savetxt(f"{sv_name}_{target_text}_dof.csv", dof)

            print(f'Leakage was detected in {sum(leakage2)} points')

            self.chi_diagrams( p, target_text, neg_log_v)

            if self.return_results:
                return None, None, chi2, p, dof, expected_values, leakage2


        elif opt == 2:
            self.file_x = self.file[f"{self.keys[self.x]}"]

            print("Starting Welch's t-test")
            sv_name = self.name.replace(".h5", "")

            result, leakage = self.start_t_test()
            target_text = str(target_byte)
            np.savetxt(f"{sv_name}_{self.keys[self.x]}_{target_text}_t-test_t.csv", result)

            print(f'Leakage was detected in {sum(leakage)} points')

            print("Starting Pearsons's chi2-test")
            self.t_diagrams(result, leakage, target_text)

            chi2, p, dof, expected_values, leakage2 = self.start_chi2_test()

            if self.generate_csv:
                np.savetxt(f"{sv_name}_{self.keys[self.x]}_{target_text}_Chi2.csv", chi2)
                np.savetxt(f"{sv_name}_{self.keys[self.x]}_{target_text}_Chi2_p.csv", p)
                np.savetxt(f"{sv_name}_{self.keys[self.x]}_{target_text}_dof.csv", dof)

            print(f'Leakage was detected in {sum(leakage2)} points')

            if self.draw_diagram:
                self.chi_diagrams(p, target_text, neg_log_v)

                plot_name = self.name.replace(".h5", "")
                plt.title(f"Comparing Welch's t-test p-values and Pearson's chi2-test p-values - {plot_name} ")
                for count, res in enumerate(result):
                    for cc, re in enumerate(res):
                        if re < 0:
                            res[cc] = -re
                    if count == 0:
                        plt.plot(res, label="|t| values", color='green')
                    else:
                        plt.plot(res, color='green')
                lneg = neg_log_v(p)
                plt.plot(lneg, label="-log10(p)", color='purple')
                plt.axhline(y=4.5, color='red')
                plt.xlabel("measurement times")
                plt.ylabel("|t|-values and -log10(p)")
                plt.legend()
                plot_name = self.name.replace(".h5", "")
                plot_name = plot_name + "_" + self.keys[self.x] + "_" + target_text + "_Welch's t-test VS Pearsons chi2-test"
                save_name = plot_name + ".png"
                plt.savefig(save_name)
                plt.clf()

            if self.return_results:
                return result, leakage, chi2, p, dof, expected_values, leakage2


        elif opt == 3:
            exit(0)

#This part is the menu what ensures that the user chooses an available option
    def menu(self, opt=-1, b=-1):
        if (opt > -1):
            return self.option_choosed(opt, b)

        print('''
        OPTIONS:
        [0] Leakage detection with Welch's T-test 
        [1] Leakage detection with Pearson's chi2-test
        [2] Run both tests - option [0] + option [1]
        [3] Exit
        ''')
        while (True):
            print("Please write the number of the desired option")
            try:
                opt = int(input())
                if (self.x > 3 or self.x < 0):
                    continue
                break
            except:
                continue
        return self.option_choosed(opt)

#This method converts MATLAB files into hdf5 file that is by the same system as ASCAD
#It is important to note that the matlab file has to be the same format as the AES_RD traces
    def mat_to_hdf5(self, input_file, output_file, keyy):
        mat = scipy.io.loadmat(input_file)
        f = h5py.File(output_file, "w")
        grp = f.create_group("measurements")
        grp.name
        '/measurements'

        arr = np.array(mat['CompressedTraces'], dtype=int)
        arr = np.rot90(arr)
        dset = grp.create_dataset("traces", data=arr)

        arr2 = np.array(mat['plaintext'], dtype=np.int32)
        arr3 = keyy
        dt = np.dtype([('plaintext', np.int32, (16,)), ('key', np.int32, (16,))])

        a = np.array([], dtype=dt)

        t = tqdm(total=len(arr2))
        for j, i in enumerate(arr2):
            t.update()
            data = np.array([(arr2[j], arr3)], dtype=dt)
            a = np.append(a, data)

        deset2 = grp.create_dataset("metadata", data=a)
        deset2.name
        '/measurements/metadata'

        f.close()
        print("done")
    # This method handles the trace groop selection within the selected trace file
    #There is also an option to declare here the target byte and the option for analysis
    #Thees ptions were made so the program could be easily used in jupyter notebook and
    # be integrated easily into other code
    def choose_tarce(self, x=-1, opt=-1, b=-1):
        if (x > -1):
            if (opt > -1 and b > -1):
                self.x = x
                self.target_byte = b
                self.traces = self.file[f"{self.keys[self.x]}/traces"]
                self.nb_traces = self.traces.shape[0]
                self.trace_len = self.traces.shape[1]
                self.plain = self.file[f"{self.keys[self.x]}/metadata"]["plaintext"][:, :]
                self.plain_len = len(self.plain[0])
                return self.menu(opt, b)
            else:
                self.x = x
                self.traces = self.file[f"{self.keys[self.x]}/traces"]
                self.nb_traces = self.traces.shape[0]
                self.trace_len = self.traces.shape[1]
                self.plain = self.file[f"{self.keys[self.x]}/metadata"]["plaintext"][:, :]
                self.plain_len = len(self.plain[0])
                return self.menu()

        while (True):

            print("Please write the number of the desired trace group")
            print(self.x)
            try:
                self.x = int(input())

                if (self.x > len(self.keys) - 1 or self.x < 0):
                    continue
                self.traces = self.file[f"{self.keys[self.x]}/traces"]
                self.nb_traces = self.traces.shape[0]
                self.trace_len = self.traces.shape[1]
                self.plain = self.file[f"{self.keys[self.x]}/metadata"]["plaintext"][:, :]
                self.plain_len = len(self.plain[0])
                break
            except Exception:
                continue
        return self.menu()


class RunningVariance():
    # initializes the mean, the m what is used to keep track of the sum of squares
    # also a counter n that counts the data points seen so far
    def __init__(self):
        self.mean = 0
        self.m = 0
        self.n = 0

    #  Update the running variance with a new data point and mean value.
    def update(self, x_n, mean):
        x_n = np.array(x_n, dtype=float)
        self.mean = mean
        self.m = (self.n * self.m + x_n ** 2) / (self.n + 1)
        self.n += 1
    #Returns the current running variance
    def __call__(self):
        return self.m - self.mean() ** 2


class RunningMean():
    # Initialize the mean and number of data points
    def __init__(self):
        self.mean = 0
        self.n = 0

    # Update the mean and number of data points
    def update(self, x_n):
        x_n = np.array(x_n, dtype=float)
        self.mean = (self.n * self.mean + x_n) / (self.n + 1)
        self.n += 1

    # Return the current mean
    def __call__(self):
        return self.mean


class RunningTtest():
    # The constructor initializes two running mean and variance objects
    # and a count of the number of traces in each sample.
    def __init__(self):
        self.mean = [RunningMean(), RunningMean()]
        self.variance = [RunningVariance(), RunningVariance()]
        self.number_of_traces = [0, 0]

    # This method updates the running mean and variance for the sample
    # indicated by the given label, using the new data point x.
    def update(self, x, label):
        self.mean[label].update(x)
        self.variance[label].update(x, self.mean[label])
        self.number_of_traces[label] += 1

    # This method calculates and returns the t-statistic based on the
    # current state of the running mean and variance objects.
    def __call__(self):
        t_stat = (self.mean[0]() - self.mean[1]()) / (
            np.sqrt(self.variance[0]() / self.number_of_traces[0] + self.variance[1]() / self.number_of_traces[1]))

        divident_df = (self.variance[0]() / self.number_of_traces[0] + self.variance[1]() / self.number_of_traces[
            1]) ** 2
        divisor_df = (self.variance[0]() / self.number_of_traces[0]) ** 2 / (self.number_of_traces[0] - 1) + (
                    self.variance[1]() / self.number_of_traces[1]) ** 2 / (self.number_of_traces[1] - 1)
        dof = divident_df / divisor_df

        2 * stats.t.cdf(t_stat, df=dof)

        return t_stat

#This part is a menu that scans user input determines the correctness of itthen executes the choosen function.
# the main options are to convert Matlab file to HDF5 and to run the leakage tests.
# It also reads the file name and checks its existence and correctness.
option = ""
while True:
    option = input("Enter option:\n[0]Submit HDF5 file for analysis \n[1] Convert MATLAB file to HDF5 file \n[2]Exit\n")
    if option == "0" or  option == "1":
        if option == "0":
            file_extension = ".h5"
            file_type_name = "HDF5"

        elif option == "1":
            file_extension = ".mat"
            file_type_name = "MATLAB"

        file_name = input(f"Enter a {file_type_name} file name (must end with {file_extension}) \n")
        if len(file_name) == 0:
           print("Name cannot be empty.")
        elif not file_name.endswith(file_extension):
           print(f"File name must end with {file_extension}. {file_name}")
        elif not os.path.isfile(file_name):
           print("File not found.")
        else:
            if option == "0":
               ana = trace_analyser()
               ana.name = file_name
               ana.file = h5py.File(file_name, "r")
               for i, key in enumerate(ana.file.keys()):
                   print(f'[{i}]', key)
                   ana.keys.append(key)
               ana.choose_tarce()
            elif option == "1":
               new_name = ""
               while True:
                   new_name = input("Enter the new file name without file extension \n")
                   if len(new_name) == 0:
                       print("Name cannot be empty.")
                   elif new_name.count(".") > 0:
                       print("Whitout file extension the name cannot contain a dot.")
                   else:
                       break

               key_arr = None
               while True:
                   input_str = input("Enter the key as a comma-separated list of integers: \n")
                   input_list = input_str.split(",")
                   if len(input_list) > 16 or len(input_list) < 16:
                       print("Invalid input: Input contains 16 integers.")
                       continue
                   try:
                       key_arr = np.array(input_list, dtype=np.int32)
                       break
                   except ValueError:
                       print("Invalid input: Input must contain only integers.")
                       continue

               ana = trace_analyser()
               ana.mat_to_hdf5(input_file=file_name, output_file=new_name + ".h5", keyy=key_arr)

    elif option == "2":
        exit(0)
    else:
        print("Invalid option.")

