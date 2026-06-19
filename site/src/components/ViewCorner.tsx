import Button from "@/components/Button";
import ContactForm from "@/components/ContactForm";
import Dialog from "@/components/Dialog";
import { useWindowScroll } from "@reactuses/core";
import { IconChevronUp, IconMail } from "@tabler/icons-react";

/** controls that float in corner of screen */
export default function ViewCorner() {
  const { y } = useWindowScroll();

  return (
    <div className="fixed right-2 bottom-2 flex flex-col gap-1 text-lg [&>button]:bg-white">
      {y > 100 && (
        <Button
          design="bubble"
          className="p-0!"
          data-tooltip="Jump to top"
          onClick={() => window.scrollTo({ top: 0 })}
        >
          <IconChevronUp />
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
            <IconMail />
          </Button>
        }
      >
        <ContactForm />
      </Dialog>
    </div>
  );
}
