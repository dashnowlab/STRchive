import type { ReactNode } from "react";
import Button from "@/components/Button";
import Popover from "@/components/Popover";
import { IconHelpCircle } from "@tabler/icons-react";

type Props = {
  children?: ReactNode;
};

export default function Help({ children }: Props) {
  if (!children) return null;
  return (
    <Popover content={children}>
      <Button>
        <IconHelpCircle className="text-gray" />
      </Button>
    </Popover>
  );
}
