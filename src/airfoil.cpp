#include "op_seq.h"

#include "cgnslib.h"

int main(int argc, char **argv) {
  // Initialise OP2
  op_init(argc, argv, 2);

  // Read CGNS grid
  int file;
  if(cg_open("./naca0012.cgns", CG_MODE_READ, &file)) {
    cg_error_exit();
  }

  cg_close(file);

  // Clean up OP2
  op_exit();
}
