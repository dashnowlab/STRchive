import { IconHelpCircle } from "@tabler/icons-react";

type Props = {
  children: string;
};

export default function Help({ children }: Props) {
  return <IconHelpCircle className="text-gray" data-tooltip={children} />;
}
