import { FaAngleUp, FaRegEnvelope } from "react-icons/fa";
import { useWindowScroll } from "@reactuses/core";
import Dialog from "@/components/Dialog";
import ContactForm from "./ContactForm";
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
        title="Contact Us"
        trigger={
          <button type="button" className="button" data-tooltip="Contact us">
            <FaRegEnvelope />
          </button>
        }
      >
        <ContactForm />
      </Dialog>
    </div>
  );
};

export default ViewCorner;
