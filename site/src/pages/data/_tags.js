import { BsStars } from "react-icons/bs";
import { FaCircleExclamation } from "react-icons/fa6";

export const tagOptions = [
  {
    value: "conflicting",
    Icon: FaCircleExclamation,
    color: "var(--error)",
  },
  {
    value: "new",
    Icon: BsStars,
    color: `var(--primary)`,
  },
];
