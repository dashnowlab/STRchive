import { LuCircleHelp } from "react-icons/lu";

type Props = {
  children: string;
};

export default function Help({ children }: Props) {
  return <LuCircleHelp className="text-gray" data-tooltip={children} />;
}
