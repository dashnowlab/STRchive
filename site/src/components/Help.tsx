import type { ReactNode } from "react";
import Popover from "@/components/Popover";
import { IconHelpCircle } from "@tabler/icons-react";

type Props = {
  children?: ReactNode;
};

export default function Help({ children }: Props) {
  if (!children) return null;
  return (
    <Popover content={children}>
      <button type="button">
        <IconHelpCircle className="text-gray" />
      </button>
    </Popover>
  );
}
