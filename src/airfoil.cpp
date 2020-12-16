#include "op_seq.h"

int main(int argc, char **argv) {
  // Initialise OP2
  op_init(argc, argv, 2);

  op_printf("Minimal OP2 Program\n");

  // Clean up OP2
  op_exit();
}
