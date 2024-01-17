
extern "C" void TOVX_write_1D_datafile(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_TOVX_write_1D_datafile
  DECLARE_CCTK_PARAMETERS
  int i;

  if (TOV_Num_TOVs > 1)
    CCTK_WARN(0, "Writing a data file for multiple stars is not (yet) "
                 "supported");
  FILE *file;
  file = fopen(TOV_save_to_datafile, "w");
  fprintf(file, "TOVSolverX data file\n");
  fprintf(file, "version 1.0\n");
  fprintf(file, "TOV_Num_Radial %d\n", (int)TOV_Num_Radial);
  fprintf(file, "\n");
  for (i=0; i < TOV_Num_Radial; i++)
  {
      fprintf(file, "%g %g %g %g %g\n",
              TOV_rbar_1d[i],
              TOV_r_1d[i],
              TOV_m_1d[i],
              TOV_phi_1d[i],
              TOV_press_1d[i]);
  }
  fclose(file);
  CCTK_WARN(0, "Simulation stopped as requested after writing data to file");
}
