import Button from "@/components/Button";
import ContactForm from "@/components/ContactForm";
import Dialog from "@/components/Dialog";

export default function ContactLink() {
  return (
    <Dialog
      title="Contact Us"
      trigger={
        <Button
          id="footer-contact"
          className="underline"
          data-tooltip="Contact us"
        >
          Contact
        </Button>
      }
    >
      <ContactForm />
    </Dialog>
  );
}
