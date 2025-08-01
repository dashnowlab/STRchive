import { FaAngleUp, FaComment } from "react-icons/fa";
import { useWindowScroll } from "@reactuses/core";
import Dialog from "@/components/Dialog";
import Feedback from "./Feedback";
import classes from "./ViewCorner.module.css";

/** controls that float in corner of screen */
const ViewCorner = () => {
  const { y } = useWindowScroll();
  return (
    <div className={classes.corner}>
      {y > 100 && (
        <button
          className="button"
          data-tooltip="Jump to top"
          onClick={() => window.scrollTo({ top: 0 })}
        >
          <FaAngleUp />
        </button>
      )}
      <Dialog
        title="Feedback"
        trigger={
          <button
            type="button"
            className="button"
            data-tooltip="Give us feedback"
          >
            <FaComment />
          </button>
        }
      >
        <Feedback />
      </Dialog>
    </div>
  );
};

export default ViewCorner;
