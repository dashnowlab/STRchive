import { FaAngleUp, FaRegEnvelope } from "react-icons/fa6";
import { useWindowScroll } from "@reactuses/core";
import Dialog from "@/components/Dialog";
import Button from "@/components/Button";
import ContactForm from "@/components/ContactForm";
import classes from "./ViewCorner.module.css";

/** controls that float in corner of screen */
const ViewCorner = () => {
  const { y } = useWindowScroll();

  return (
    <div className={classes.corner}>
      {y > 100 && (
        <Button
          design="bubble"
          data-tooltip="Jump to top"
          onClick={() => window.scrollTo({ top: 0 })}
        >
          <FaAngleUp />
        </Button>
      )}
      <Dialog
        title="Contact Us"
        trigger={
          <Button design="bubble" id="corner-contact" data-tooltip="Contact us">
            <FaRegEnvelope />
          </Button>
        }
      >
        <ContactForm />
      </Dialog>
    </div>
  );
};

export default ViewCorner;
