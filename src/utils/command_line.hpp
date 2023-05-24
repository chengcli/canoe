#ifndef SRC_UTILS_COMMAND_LINE_HPP_
#define SRC_UTILS_COMMAND_LINE_HPP_

struct CommandLine {
  static char *input_filename;
  static char *restart_filename;
  static char *prundir;
  static int res_flag;
  static int narg_flag;
  static int iarg_flag;
  static int mesh_flag;
  static int wtlim;
  static int argc;
  static char **argv;

  CommandLine(int argc, char **argv);
  ~CommandLine() {}
};

#endif  // SRC_UTILS_COMMAND_LINE_HPP_
