import { FaRegCircleQuestion } from "react-icons/fa6";

const Help = ({ children }) => (
  <FaRegCircleQuestion className="text-gray" data-tooltip={children} />
);

export default Help;
