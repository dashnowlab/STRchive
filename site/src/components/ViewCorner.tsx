import { LuChevronUp, LuMail } from "react-icons/lu";
import Button from "@/components/Button";
import ContactForm from "@/components/ContactForm";
import Dialog from "@/components/Dialog";
import { useWindowScroll } from "@reactuses/core";

/** controls that float in corner of screen */
export default function ViewCorner() {
  const { y } = useWindowScroll();

  return (
    <div className="fixed right-2.5 bottom-2.5 flex flex-col gap-1 text-lg [&>button]:bg-white">
      {y > 100 && (
        <Button
          design="bubble"
          className="p-0!"
          data-tooltip="Jump to top"
          onClick={() => window.scrollTo({ top: 0 })}
        >
          <LuChevronUp />
        </Button>
      )}
      <Dialog
        title="Contact Us"
        trigger={
          <Button
            design="bubble"
            className="p-0!"
            id="corner-contact"
            data-tooltip="Contact us"
          >
            <LuMail />
          </Button>
        }
      >
        <ContactForm />
      </Dialog>
    </div>
  );
}
