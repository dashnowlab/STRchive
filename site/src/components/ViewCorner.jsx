import { FaAngleUp, FaRegEnvelope } from "react-icons/fa6";
import { useWindowScroll } from "@reactuses/core";
import Button from "@/components/Button";
import ContactForm from "@/components/ContactForm";
import Dialog from "@/components/Dialog";

/** controls that float in corner of screen */
const ViewCorner = () => {
  const { y } = useWindowScroll();

  return (
    <div className="fixed bottom-2.5 right-2.5 flex flex-col gap-[5px] text-[1.1rem] [&>button]:bg-white">
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
