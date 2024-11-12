import { BsStars } from "react-icons/bs";
import { FaCircleExclamation } from "react-icons/fa6";
import { newThreshold } from "./_derived";

/** top-level tag types */
export const tagOptions = [
  {
    value: "conflicting",
    Icon: FaCircleExclamation,
    color: "var(--secondary)",
    tooltip: "Conflicting evidence",
  },
  {
    value: "new",
    Icon: BsStars,
    color: `var(--primary)`,
    tooltip: `Relatively new (less than ~${newThreshold} years old)`,
  },
];
