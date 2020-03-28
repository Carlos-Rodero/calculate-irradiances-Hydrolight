from calculate_irradiances import ProcessRadFile


if __name__ == "__main__":

    prf = ProcessRadFile()
    prf.calc_irradiances()
    prf.plot_irradiances()
