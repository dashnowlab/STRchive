import ContactForm from "@/components/ContactForm";
import Dialog from "@/components/Dialog";

export default function ContactLink() {
  return (
    <Dialog
      title="Contact Us"
      trigger={
        <button id="footer-contact" className="underline">
          Contact
        </button>
      }
    >
      <ContactForm />
    </Dialog>
  );
}
